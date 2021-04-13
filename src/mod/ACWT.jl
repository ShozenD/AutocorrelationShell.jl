module ACWT
export
    convf,
    iconv,
    echant,
    acfilter,
    leftchild,
    rightchild,
    parentnode,
    subtree

using SpecialFunctions
using Wavelets, LinearAlgebra

## AC1D & AC2D
"""Convolution filter"""
# function convf(a,b,sig)
#     rst = zeros(size(sig))
#     for i in eachindex(sig)
#         tmp =  sum(b[j] * sig[i-j+1] for j in 1:min(i, length(b)))
#         tmp -= sum(a[j] * rst[i-j+1] for j in 1:min(i, length(a)))
#         rst[i] = tmp / a[1]
#     end
#     return rst
# end

# function iconv(f::Vector{T},x::Vector{T}) where T <: Number
#    n = length(x)
#    p = length(f)
#    if p <= n
#       xpadded = vcat(x[(n+1-p):n], x)
#    else
#       z = zeros(p)
#       for i=1:p
#           imod = 1 + rem(p*n -p + i-1,n)
#           z[i] = x[imod]
#       end
#       xpadded = vcat(z,x)
#    end
#    ypadded = convf(1, f, xpadded)
#    y = ypadded[(p+1):(n+p)]
#    return y
# end

"""Sample the range `(b+1):n` in intervals of `2^d`"""
# function echant(n::Integer, d::Integer, b::Integer)
#     return (b + 1):(2^d):n
# end

"""Computes ac coefficients of given filter"""
function autocorr(f::OrthoFilter)
    H = WT.qmf(f)
    l = length(H)
    result = zeros(l - 1)
    for k in 1:(l - 1)
        for i in 1:(l - k)
            @inbounds result[k] += H[i] * H[i + k]
        end
        result[k] *= 2
    end
    return result
end

"""Computes the response of signal to the autocorrelation filter"""
# function acfilter(x::Vector{T}, f::Vector{T}) where T <: Number
#     n = length(x)
#     p = length(f)
#     tran2 = p - 1
#     tran1 = tran2 ÷ 2
# 
#     d = circshift(iconv(f, circshift(x, tran1)), -tran2)
#     return d[1:n]
# end

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

function makereverseqmfpair(f::OrthoFilter)
    pmf, qmf = makeqmfpair(f)
    return reverse(pmf), reverse(qmf)
end

"""For acwpt array methods"""
leftchild(i::Integer) = i<<1;
rightchild(i::Integer) = i<<1+1;
children(i::Integer) = (leftchild(i),rightchild(i));
parentnode(i::Integer) = i>>1;

"""
    subtree(i, M)

Given the index `i` of a node in a binary tree, find all the index numbers of its children.
`M` is the maximum possible index number for that tree (the bottom right node).
"""
function subtree(i::Integer, M::Integer)
    nodes = Vector{Integer}(undef,0)
    subtree!(nodes,i,M)
    return nodes
end

function subtree!(nodes::Vector{Integer}, i::Integer, M::Integer)
    if i<<1 <= M
        append!(nodes, children(i))
        subtree!(nodes, leftchild(i), M)
        subtree!(nodes, rightchild(i), M)
    end
end

end # End module
