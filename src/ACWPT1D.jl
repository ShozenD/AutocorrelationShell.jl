## Define base node object
mutable struct BinaryNode{T}
    data::T
    parent::BinaryNode{T}
    left::BinaryNode{T} # pointers to children
    right::BinaryNode{T}

    # for graphing purposes
    depth::Int # tree depth

    # Root constructor
    BinaryNode{T}(data) where T = new{T}(data)
    # Child node constructor
    BinaryNode{T}(data, parent::BinaryNode{T}) where T = new{T}(data, parent)
end

BinaryNode(data) = BinaryNode{typeof(data)}(data)

## base node methods
function initializeBinaryNode(data) # Rootnode initialize method
    node = BinaryNode(data)
    node.depth = 0
    return node
end

function initializeBinaryNode(data, parent::BinaryNode, depth::Int) # childnode initialize method
    node = typeof(parent)(data, parent)
    node.depth = depth
    return node
end

function leftchild(data, parent::BinaryNode, depth::Int)
    !isdefined(parent, :left) || error("left child is already assigned")
    node = initializeBinaryNode(data, parent, depth)
    parent.left = node
end

function rightchild(data, parent::BinaryNode, depth::Int)
    !isdefined(parent, :right) || error("right child is already assigned")
    node = initializeBinaryNode(data, parent, depth)
    parent.right = node
end

## Using the AbstractTrees API
function AbstractTrees.children(node::BinaryNode)
    if isdefined(node, :left)
        if isdefined(node, :right)
            return (node.left, node.right)
        end
        return (node.left,)
    end
    isdefined(node, :right) && return (node.right,)
    return ()
end

# Things that make the printing easier
AbstractTrees.printnode(io::IO, node::BinaryNode) = print(io, node.data)

## Optional enhancements
# These next two definitions allow inference of the item type in iteration.
# (They are not sufficient to solve all internal inference issues, however.)
Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = BinaryNode{T}
Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = Base.HasEltype()

## Wavelet Packets functions
function acwpt(x::Vector{T}, node::BinaryNode,
    P::Vector{T}, Q::Vector{T}, d::Integer) where T<:Real

    n = length(x)
    J = dyadlength(x) # maximum tree depth
    left, right = zeros(n), zeros(n)

    if d < J
        @inbounds begin
            for b = 0:(2^d-1) # depth starts from 2
                s = x[echant(n, d, b)]
                h = ac_filter(s, Q)
                l = ac_filter(s, P)
                left[echant(n, d, b)] = l
                right[echant(n, d, b)] = h
            end
        end

        # left
        leftchild(left, node, d+1)
        acwpt(left, node.left, P, Q, d+1)

        # right
        rightchild(right, node, d+1)
        acwpt(right, node.right, P, Q, d+1)
    end
end

"""
    acwpt(x, P, Q)

Compute the autocorrelation wavelet packet transform for a given signal. Returns a binary tree object.

# Arguments
- `x::Vector{<:Real}`: 1 dimensional signal.
- `P::Vector{<:Real}`: Low autocorrelation shell filter.
- `Q::Vector{<:Real}`: High autocorrelation shell filter.
"""
function acwpt(x::Vector{T}, P::Vector{T}, Q::Vector{T}) where T<:Real
    root = initializeBinaryNode(x) # original signal
    acwpt(x, root, P, Q, 0)
    return root
end

"""
    iacwpt(tree)

Reconstructs the signal using the autocorrelation wavelet packet trasform bases.

# Arguments
- `tree::BinaryNode`: The root node of the wavelet packet tree (binary tree)
"""
function iacwpt(tree::BinaryNode)
    if isdefined(tree, :left)
        left = iacwpt(tree.left)
    end

    if isdefined(tree, :right)
        right = iacwpt(tree.right)
    end

    if !isdefined(tree, :left) & !isdefined(tree, :right)
        return tree.data
    end

    return (left + right)/sqrt(2)
end

## Post Order Best Basis
function acwptPostOrderBestBasis!(node::BinaryNode; direction::AbstractString="right", et::Wavelets.Entropy=NormEntropy())
    if isdefined(node, :left) & isdefined(node, :right)

        el = acwptPostOrderBestBasis!(node.left, direction="left")
        er = acwptPostOrderBestBasis!(node.right)

        if (el + er)/2 < wentropy(node.data, et)
            return (el + er)/2

        elseif isdefined(node, :parent) # handling cases where the first decomposition is unuseful
            data = node.data
            # new_pruned_node = typeof(node)(data, node) # create node with no children
            new_pruned_node = initializeBinaryNode(data, node, node.depth)
            if direction == "right"
                node.parent.right = new_pruned_node
            else
                node.parent.left = new_pruned_node
            end
            return wentropy(node.data, et)
        end
    else
        return wentropy(node.data, et) # return entropy value
    end
end

"""
    acwptPostOrderBestBasis(tree::BinaryNode; et::Wavelets.Entropy=NormEntropy())

Returns the best basis tree found using post order traversal. This is a democratic approach to finding the best basis tree.

# Arguments
-`tree::BinaryNode`: The entry point of the wavelet packet decomposition tree.
-`et::Wavelets.Entropy`: The type of cost function used for evaluation.
"""
function acwptPostOrderBestBasis(tree::BinaryNode; et::Wavelets.Entropy=NormEntropy())
    best_tree = deepcopy(tree)
    acwptPostOrderBestBasis!(best_tree, et=et)
    return best_tree
end

function acwptPreOrderBestBasis!(node::BinaryNode; direction::AbstractString="right", et::Wavelets.Entropy=NormEntropy())
    if isdefined(node, :left) & isdefined(node, :right)
        e0 = wentropy(node.data, et)
        e10 = wentropy(node.left.data, et)
        e11 = wentropy(node.right.data, et)
        if (e10 + e11)/2 < e0
            acwptPreOrderBestBasis!(node.left, direction="left", et=et)
            acwptPreOrderBestBasis!(node.right, et=et)
        else
            data = node.data
            new_pruned_node = typeof(node)(data, node) # create node with no children
            if isdefined(node, :parent) # handling cases where the first decomposition is unuseful
                if direction == "right"
                    node.parent.right = new_pruned_node
                else
                    node.parent.left = new_pruned_node
                end
            else print("No good tree available")
            end
        end
    end
end


"""
    acwptPreOrderBestBasis(node::BinaryNode; et::Wavelets.Entrop=NormEntropy())

Returns the best basis tree found using pre order traversal. This is a greedy approach to finding the best basis tree.

# Arguments
-`tree::BinaryNode`: The entry point of the wavelet packet decomposition tree.
-`et::Wavelets.Entropy`: The type of cost function used for evaluation.
"""
function acwptPreOrderBestBasis(tree::BinaryNode; et::Wavelets.Entropy=NormEntropy())
    best_tree = deepcopy(tree)
    acwptPreOrderBestBasis!(best_tree)
    return best_tree
end

"""
    nodes_to_bitmap(x::BinaryNode)

Translates the binary tree datastructure to a binary bitmap in preparation for visualization.

# Arguments
-`x::BinaryNode`: The entry point of the best basis tree.
"""
function tree_to_bitmap(x::BinaryNode)
    nrow, ncol = length(x.data), dyadlength(x.data) + 1
    arr, hash = zeros(nrow, ncol), zeros(ncol - 1)

    @inbounds begin
        for n in collect(PreOrderDFS(x))
            d = n.depth
            l = length(n.data)/2^d
            _start = Int(hash[d]+1)
            _end = Int(_start + l - 1)
            if !isdefined(n, :left)
                hash[d+1:end] .= _start + l/2 - 1
                arr[_start:_end, d+1] .= 1
            end
            if !isdefined(n, :right)
                hash[d+1:end] .= _end
                arr[_start:_end, d+1] .= 1
            end
            hash[d] += l
        end
    end

    return arr
end

"""
    plot_tfbdry(x::BinaryNode)

Plots the selected nodes of the best basis tree

# Arguments
-`x::BinaryNode`: The entry point of the best basis tree
"""
function plot_tfbdry(x::BinaryNode)
    bitmap = tree_to_bitmap(x)
    nrow, ncol = size(bitmap)
    heatmap(
        transpose(bitmap),
        yflip=true,
        legend=false,
        yticks=([1:1:(ncol);], map(string, 0:1:(ncol-1)))
    )
    hline!([1.5:1:(ncol-0.5);], color=:red)
    @inbounds begin
        for i in 1:ncol
            for j in 1:2^(i-1)
                vpos = (nrow/2^i)*(2*j-1) + 0.5
                plot!(vpos*ones(ncol-i+1), (i+0.5):(ncol+0.5), color=:red)
            end
        end
    end
    current()
end
