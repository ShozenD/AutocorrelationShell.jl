module ACTransforms
export
    acwt,
    iacwt,
    acwpt,
    iacwpt
using ..ACWT
using LinearAlgebra, Wavelets

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
	iacwt

The inverse of `acwt`

**See also:** `acwt`
"""
function iacwt end

"""
    iacwpt

The inverse of `acwpt`

**See also:** `acwpt`
"""

function iacwt!(x::AbstractArray{<:Number,2})
    n,m = size(x)
    for i = 2:m
        @inbounds x[:,1] = (x[:,1] + x[:,i])/sqrt(2)
    end
end

function iacwt(x::AbstractArray{<:Number,2})
    y = deepcopy(x)
    iacwt!(y)
    return y[:,1]
end

function acwt_step(v::AbstractVector{T}, j::Integer, h::Array{T,1}, g::Array{T,1}) where T <: Number
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
        @inbounds wp[:,1], wp[:,L+2-j] = acwt_step(wp[:,1],j,Qmf,Pmf)
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

function acwpt_step(W::AbstractArray{T,2}, i::Integer, d::Integer, Qmf::Vector{T}, Pmf::Vector{T}) where T <: Number
    n,m = size(W)
    if i<<1+1 <= m
      W[:,i<<1], W[:,i<<1+1] = acwt_step(W[:,i],d,Qmf,Pmf)
      acwpt_step(W, leftchild(i), d+1, Qmf, Pmf) # left
      acwpt_step(W, rightchild(i), d+1, Qmf, Pmf) # right
    end
  end
  
  function acwpt(x::AbstractVector{T}, wt::OrthoFilter, L::Integer=maxtransformlevels(x)) where T<:Number
    W = Array{Float64,2}(undef,length(x),2<<L-1)
    W[:,1] = x
    Pmf, Qmf = ACWT.makereverseqmfpair(wt)
    acwpt_step(W,1,1,Qmf,Pmf)
    return W
  end

function iacwpt(xw::AbstractArray{<:Number,2}, tree::BitVector, i::Integer=1)

    @assert i <= size(xw, 2)
    @assert isvalidtree(xw[:,1], tree)
    n₀ = length(tree)
    if i > n₀ || tree[i] == false       # leaf node
        return xw[:,i]
    end

    v₀ = iacwpt(xw,tree,leftchild(i))
    v₁ = iacwpt(xw,tree,rightchild(i))

    return (v₀ + v₁) / √2
end

end # End module
