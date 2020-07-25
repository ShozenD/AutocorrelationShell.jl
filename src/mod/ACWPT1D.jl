module ACWPT1D
export
    acwpt

using AbstractTrees
using ..AC1D

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

end # module
