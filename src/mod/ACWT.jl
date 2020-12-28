module ACWT
export
    convf,
    iconv,
    echant,
    acfilter,
    AcwptNode,
    initAcwptNode,
    leftchild,
    rightchild

using SpecialFunctions
using AbstractTrees, Wavelets, LinearAlgebra

## AC1D & AC2D
"""Convolution filter"""
function convf(a,b,sig)
    rst = zeros(size(sig))
    for i in eachindex(sig)
        tmp =  sum(b[j] * sig[i-j+1] for j in 1:min(i, length(b)))
        tmp -= sum(a[j] * rst[i-j+1] for j in 1:min(i, length(a)))
        rst[i] = tmp / a[1]
    end
    return rst
end

function iconv(f::Vector{T},x::Vector{T}) where T <: Number
   n = length(x)
   p = length(f)
   if p <= n
      xpadded = vcat(x[(n+1-p):n], x)
   else
      z = zeros(p)
      for i=1:p
          imod = 1 + rem(p*n -p + i-1,n)
          z[i] = x[imod]
      end
      xpadded = vcat(z,x)
   end
   ypadded = convf(1, f, xpadded)
   y = ypadded[(p+1):(n+p)]
   return y
end

"""Sample the range `(b+1):n` in intervals of `2^d`"""
function echant(n::Integer, d::Integer, b::Integer)
    return (b + 1):(2^d):n
end

"""Computes ac coefficients of given filter"""
function autocorr(f::OrthoFilter)
    H = WT.qmf(f)
    l = length(H)
    result = zeros(l - 1)
    @inbounds begin
        for k in 1:(l - 1)
            for i in 1:(l - k)
                result[k] += H[i] * H[i + k]
            end
            result[k] *= 2
        end
    end
    return result
end

"""Computes the response of signal to the autocorrelation filter"""
function acfilter(x::Vector{T}, f::Vector{T}) where T <: Number
    n = length(x)
    p = length(f)
    tran2 = p - 1
    tran1 = tran2 ÷ 2

    d = circshift(iconv(f, circshift(x, tran1)), -tran2)
    return d[1:n]
end

"""Computes low ac filter"""
function pfilter(f::OrthoFilter)
    a = autocorr(f)
    c1 = 1 / sqrt(2)
    c2 = c1 / 2
    b = c2 * a
    return vcat(b[end:-1:1], c1, b)
end

"""Computes the high AC filter"""
function qfilter(f::OrthoFilter)
    a = autocorr(f)
    c1 = 1 / sqrt(2)
    c2 = c1 / 2
    b = -c2 * a
    return vcat(b[end:-1:1], c1, b)
end

function makeqmfpair(f::OrthoFilter)
    pmfilter, qmfilter = pfilter(f), qfilter(f)
    return pmfilter, qmfilter
end

## ACWPT
"""
The `AcwptNode` type is a composite type to records the result of the
ac wavelet packet transform and its meta data.
"""
abstract type Node end
abstract type BinaryNode <: Node end
mutable struct AcwptNode{T} <: BinaryNode
    data::T
    parent::AcwptNode{T}
    left::AcwptNode{T} # pointers to children
    right::AcwptNode{T}

    # for graphing purposes
    depth::Int # tree depth

    # Root constructor
    AcwptNode{T}(data) where T = new{T}(data)
    # Child node constructor
    AcwptNode{T}(data, parent::AcwptNode{T}) where T = new{T}(data, parent)
end

"""Initializes root AcwptNode"""
function AcwptNode(x)
    node = AcwptNode{typeof(x)}(x)
    node.depth = 0
    return node
end

"""Initialize child AcwptNode"""
function AcwptNode(x, parent::AcwptNode, depth::Integer=parent.depth+1)
    node = typeof(parent)(x,parent)
    node.depth = depth
    return node
end

"""Creates left child node"""
function leftchild(x, parent::AcwptNode, depth::Integer=parent.depth+1)
    !isdefined(parent, :left) || error("left child is already assigned")
    parent.left = AcwptNode(x, parent, depth)
end

"""Creates right child node"""
function rightchild(data, parent::AcwptNode, depth::Integer=parent.depth+1)
    !isdefined(parent, :right) || error("right child is already assigned")
    parent.right = AcwptNode(data, parent, depth)
end

"""Using the AbstractTree API"""
function AbstractTrees.children(node::AcwptNode)
    if isdefined(node, :left)
        if isdefined(node, :right)
            return (node.left, node.right)
        end
        return (node.left,)
    end
    isdefined(node, :right) && return (node.right,)
    return ()
end

"""Things that make the printing easier"""
AbstractTrees.printnode(io::IO, node::AcwptNode) = print(io, node.data)

"""
**Optional enhancements**
These next two definitions allow inference of the item type in iteration.
(They are not sufficient to solve all internal inference issues, however.)
"""
Base.eltype(::Type{<:TreeIterator{AcwptNode{T}}}) where T = AcwptNode{T}
Base.IteratorEltype(::Type{<:TreeIterator{AcwptNode{T}}}) where T = Base.HasEltype()

end # End module
