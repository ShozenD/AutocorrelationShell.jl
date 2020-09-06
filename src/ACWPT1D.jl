## Define base node object
mutable struct BinaryNode{T}
    data::T
    parent::BinaryNode{T}
    left::BinaryNode{T}
    right::BinaryNode{T}

    # Root constructor
    BinaryNode{T}(data) where T = new{T}(data)
    # Child node constructor
    BinaryNode{T}(data, parent::BinaryNode{T}) where T = new{T}(data, parent)
end

BinaryNode(data) = BinaryNode{typeof(data)}(data)

## base node methods
function leftchild(data, parent::BinaryNode)
    !isdefined(parent, :left) || error("left child is already assigned")
    node = typeof(parent)(data, parent)
    parent.left = node
end

function rightchild(data, parent::BinaryNode)
    !isdefined(parent, :right) || error("right child is already assigned")
    node = typeof(parent)(data, parent)
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
function binary_ac1d(x::Vector{T}, node::BinaryNode,
    L::Integer, P::Vector{T}, Q::Vector{T}, depth::Integer) where T<:Real

    max_depth = L+2 # maximum tree depth

    if depth <= max_depth
        decomp = fwt_ac(x, L, P, Q)

        # left
        left = decomp[:,1]
        leftchild(left, node)
        binary_ac1d(left, node.left, L, P, Q, depth+1)

        # right
        right = decomp[:,2]
        rightchild(right, node)
        binary_ac1d(right, node.right, L, P, Q, depth+1)
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

    max_depth = dyadlength(x) + 1
    L = max_depth - 2

    root = BinaryNode(x) # original signal
    binary_ac1d(x, root, L, P, Q, 2)

    return root
end

"""
    aciwpt(tree)

Reconstructs the signal using the autocorrelation wavelet packet trasform bases.

# Arguments
- `tree::BinaryNode`: The root node of the wavelet packet tree (binary tree)
"""
function aciwpt(tree::BinaryNode)
    if isdefined(tree, :left)
        left = aciwpt(tree.left)
    end

    if isdefined(tree, :right)
        right = aciwpt(tree.right)
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
        else isdefined(node, :parent) # handling cases where the first decomposition is unuseful
            data = node.data
            new_pruned_node = typeof(node)(data, node) # create node with no children
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
