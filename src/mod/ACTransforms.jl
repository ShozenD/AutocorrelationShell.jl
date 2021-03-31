module ACTransforms
export
    acwt,
    iacwt,
    hacwt,
    vacwt,
    acwpt,
    iacwpt
using ..ACWT, ..ACUtil
using AbstractTrees, LinearAlgebra, Wavelets

"""
	acwt(x, wt[, L=maxtransformlevels(x)])

Perform a forward ac wavelet transform of the array `x`.
This method works for the 2D case as well.
The wavelet type `wt` determines the transform type.
Refer to Wavelet.jl for a list of available methods.

# Examples
```julia
acwt(x, wavelet(WT.db4))
```

**See also:** `iacwt`
"""
function acwt end

"""
    acwpt(x, wt[, L=maxtransformlevels(x)])

Perform a ac wavelet packet transform of the array `x`.
The wavelet type `wt` determines the transform type.

# Examples
```julia
acwpt(x, wavelet(WT.db4))
```

**See also:** `iacwpt`
"""
function acwpt end

"""
    bestbasistree(x, et)

Returns the best basis tree for the ac wavelet packet transform
using a specified cost function. *Default*: L1 norm (Sparcity)

# Examples
```julia
bestbasistree(acwpt(x, wavelet(WT.db4)), NormEntropy())
```

**See also:** `acwpt`, `iacwpt`, `bestbasistree!`
"""
function bestbasistree end

"""
	iacwt

The inverse of `acwt`

**See also:** `acwt`
"""
function iacwt end

"""
    iacwt!

Same as `iacwt`, but without array allocation

**See also:** `iacwt`
"""
function iacwt! end

"""
    iacwpt

The inverse of `acwpt`

**See also:** `acwpt`
"""

function iacwt!(x::AbstractArray{<:Number,2})
    n,m = size(x)
    @inbounds begin
        for i = 2:m
            x[:,1] = (x[:,1] + x[:,i])/sqrt(2)
        end
    end
end

function iacwt(x::AbstractArray{<:Number,2})
    y = deepcopy(x)
    iacwt!(y)
    return y[:,1]
end

# function acwt(x::AbstractArray{<:Number,1}, wt::OrthoFilter, L::Integer=maxtransformlevels(x))
#     # Setup
#     n = length(x)
#     Pmf, Qmf = ACWT.makeqmfpair(wt)
#     wp = zeros(n,L+1)
#     wp[:,1] = x
#     @inbounds begin
#         for d=0:(L-1)
#             for b=0:(2^d-1)
#                 s = wp[echant(n,d,b),1]
#                 high = acfilter(s,Qmf)
#                 low = acfilter(s,Pmf)
#                 wp[echant(n,d,b),L+1-d] = high
#                 wp[echant(n,d,b),1] = low
#             end
#         end
#     end
#     return wp
# end

function acwt(v::AbstractVector{T}, j::Integer, h::Array{S,1},
    g::Array{S,1}) where {T <: Number, S <: Number}
    N = length(v)
    L = length(h)
    v1 = zeros(T, N)
    w1 = zeros(T, N)

    @inbounds begin
        for i in 1:N
            t = i
            i = mod1(i + (L÷2+1) * 2^(j-1),N) # Need to shift by half the filter size because of periodicity assumption 
            for n in 1:L
                t += 2^(j-1)
                t = mod1(t, N)
                w1[i] += h[n] * v[t]
                v1[i] += g[n] * v[t]
            end
        end
    end
    return v1, w1 
end

function acwt(x::AbstractVector{<:Number}, wt::OrthoFilter, L::Integer=maxtransformlevels(x))

    @assert L <= maxtransformlevels(x) || throw(ArgumentError("Too many transform levels (length(x) < 2^L"))
    @assert L >= 1 || throw(ArgumentError("L must be >= 1"))

    # Setup
	n = length(x)
    Pmf, Qmf = ACWT.makereverseqmfpair(wt)
	wp = zeros(n,L+1)
	wp[:,1] = x

    for j in 1:L
        @inbounds wp[:,1], wp[:,L+2-j] = acwt(wp[:,1],j,Qmf,Pmf)
    end

	return wp
end

"""Computes the column-wise forward acwat coeficients for 2D signals."""
function hacwt(x::AbstractArray{<:Number,2},
                         wt::OrthoFilter,
                         L::Integer=maxtransformlevels(x[:,1]))
    nrow, ncol = size(x)
    W = Array{Float64,3}(undef,nrow,L+1,ncol)
    for i in 1:ncol
        @inbounds W[:,:,i] = acwt(x[:,i],wt,L)
    end
    return W
end

"""Computes the row-wise 1D acwt coeficients for 2D signals."""
function vacwt(x::AbstractArray{<:Number,2},
                      wt::OrthoFilter,
                      L::Integer=maxtransformlevels(x[1,:]))
    nrow, ncol = size(x)
    W = Array{Number,3}(undef,ncol,L+1,nrow)
    for i in 1:nrow
        W[:,:,i] = acwt(x[i,:],wt,L)
    end
    return W
end

function acwt(x::AbstractArray{<:Number,2}, wt::OrthoFilter,
              Lrow::Integer=maxtransformlevels(x[1,:]),
              Lcol::Integer=maxtransformlevels(x[:,1]))
    nrow, ncol = size(x)
    W3d = hacwt(x,wt,Lcol)
    W4d = Array{Number,4}(undef,Lcol+1,ncol,Lrow+1,nrow)
    for i in 1:Lcol+1
        @inbounds W4d[i,:,:,:] = vacwt(W3d[:,i,:],wt,Lrow)
    end
    W4d = permutedims(W4d, [4,2,3,1])
    return W4d
end

# Inverse Transforms
function iacwt(x::AbstractArray{<:Number,4})
    nrow, ncol, Lrow, Lcol = size(x)
    W4d = permutedims(x,[4,2,3,1])
    W3d = Array{Number,3}(undef, nrow, Lcol, ncol)
    for i in 1:Lcol
        for j in 1:nrow
            @inbounds W3d[j,i,:] = iacwt(W4d[i,:,:,j])
        end
    end
    y = Array{Number,2}(undef, nrow, ncol)
    for i in 1:ncol
        @inbounds y[:,i] = iacwt(W3d[:,:,i])
    end
    return y
end

## Wavelet Packets functions
# Tree method
function acwpt(x::Vector{T}, node::AcwptNode,
               wt::OrthoFilter,
               L::Integer=maxtransformlevels(x),
               d::Integer=0) where T<:Number

    n = length(x)
    Pmf, Qmf = ACWT.makeqmfpair(wt)
    left, right = zeros(n), zeros(n)

    if d < L
        @inbounds begin
            for b = 0:(2^d-1) # depth starts from 2
                s = x[echant(n,d,b)]
                high = acfilter(s,Qmf)
                low = acfilter(s,Pmf)
                left[echant(n,d,b)] = low
                right[echant(n,d,b)] = high
            end
        end

        # left
        leftchild(left,node,d+1)
        acwpt(left,node.left,wt,L,d+1)

        # right
        rightchild(right,node,d+1)
        acwpt(right,node.right,wt,L,d+1)
    end
end

# Array method
# function acwpt(W::Array{<:Number,2}, i::Integer, d::Integer,
#                Pmf::Vector{<:Number}, Qmf::Vector{<:Number})
#     n,m = size(W)
#     if i<<1+1 <= m
#         @inbounds begin
#             for b = 0:(2^d-1) # depth starts from 2
#                 s = W[echant(n,d,b),i]
#                 W[echant(n,d,b),i<<1] = acfilter(s,Pmf)
#                 W[echant(n,d,b),i<<1+1] = acfilter(s,Qmf)
#             end
#         end
#         acwpt(W,i<<1,d+1,Pmf,Qmf) # left
#         acwpt(W,i<<1+1,d+1,Pmf,Qmf) # right
#     end
# end

function acwpt(W::Array{<:Number,2}, i::Integer, d::Integer,
               Qmf::Vector{<:Number}, Pmf::Vector{<:Number})
    n,m = size(W)
    if i<<1+1 <= m
        W[:,i<<1], W[:,i<<1+1] = acwt(W[:,i],d,Qmf,Pmf)
        acwpt(W,i<<1,d+1,Qmf,Pmf) # left
        acwpt(W,i<<1+1,d+1,Qmf,Pmf) # right
    end
end

# Combined
function acwpt(x::Vector{T}, wt::OrthoFilter,
               L::Integer=maxtransformlevels(x), method::Symbol=:array) where T<:Number
    if method == :tree
        root = AcwptNode(x)
        acwpt(x,root,wt,L)
        return root
    elseif method == :array
        W = Array{Float64,2}(undef,length(x),2^(L+1)-1)
        W[:,1] = x
        Pmf, Qmf = ACWT.makereverseqmfpair(wt)
        acwpt(W,1,1,Qmf,Pmf)
        return W
    else
        throw(ArgumentError("unkown method"))
    end
end

# Inverse ac wavelet packets
function iacwpt(x::AcwptNode)
    if isdefined(x,:left)
        left = iacwpt(x.left)
    end
    if isdefined(x,:right)
        right = iacwpt(x.right)
    end
    if !isdefined(x,:left) & !isdefined(x,:right)
        return x.data
    end
    return (left+right)/sqrt(2)
end

function iacwpt(x::Array{<:Number,2}, bt::BitArray, i::Integer=1)
    M = length(bt)
    if rightchild(i) > M # Reached maximum level
        return x[:,i]
    end
    if bt[leftchild(i)] != 0
        left = iacwpt(x,bt,leftchild(i))
    end
    if bt[rightchild(i)] != 0
        right = iacwpt(x,bt,rightchild(i))
    end
    if bt[i]==1 && (bt[leftchild(i)]==0 && bt[rightchild(i)]==0) # Is leafnode
        return x[:,i]
    end
    return (left+right)/sqrt(2)
end

# TODO: Implement iacwpt for array method without inplace behavior
function iacwpt!(x::Array{<:Number,2})
    n,m = size(x)
    @inbounds begin
        for i in (m-1):-2:2
            x[:,convert(Int,i/2)] = (x[:,i] + x[:,i+1])/sqrt(2)
        end
    end
end

# best basis algorithm
function bestbasistree!(node::AcwptNode, et::Wavelets.Entropy=NormEntropy(), dir::Symbol=:right)
    if isdefined(node, :left) & isdefined(node, :right) # Is not leaf node

        el = bestbasistree!(node.left,et,:left)
        er = bestbasistree!(node.right,et,:right)

        if (el+er)/2 < Wavelets.Threshold.coefentropy(node.data,et)
            return (el+er)/2
        elseif !isdefined(node,:parent) # No interesting decomposition found
            node = AcwptNode(node.data)
        else
            if dir === :right
                node.parent.right = AcwptNode(node.data,node.parent)
            else
                node.parent.left = AcwptNode(node.data,node.parent)
            end
            return Wavelets.Threshold.coefentropy(node.data,et)
        end
    else
        return Wavelets.Threshold.coefentropy(node.data,et) # Leaf node
    end
end

function Wavelets.Threshold.bestbasistree(x::AcwptNode,et::Wavelets.Entropy=NormEntropy())
    y = deepcopy(x)
    bestbasistree!(y,et)
    return y
end

function Wavelets.Threshold.bestbasistree(W::Array{<:Number,2}, et::Wavelets.Entropy=NormEntropy())
    _, M = size(W)
    costvec = findcost(W,et)
    besttree = trues(M)

    for idx in (M>>1):-1:1
        if costvec[idx] <= (costvec[idx<<1] + costvec[(idx<<1)+1])/2
            besttree[subtree(idx,M)] .= false
        else
            costvec[idx] = (costvec[idx<<1] + costvec[(idx<<1)+1])/2
        end
    end
    return besttree
end


end # End module
