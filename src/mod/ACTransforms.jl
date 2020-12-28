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
Returns a binary x object.

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

function iacwt!(x::AbstractArray{Float64,2})
    n,m = size(x)
    @inbounds begin
        for i = 2:m
            x[:,1] = (x[:,1] + x[:,i])/sqrt(2)
        end
    end
end

function iacwt(x::AbstractArray{Float64,2})
    y = deepcopy(x)
    iacwt!(y)
    return y[:,1]
end

function acwt(x::AbstractArray{Float64,1}, wt::OrthoFilter, L::Integer=maxtransformlevels(x))
    # Setup
	n = length(x)
    Pmf, Qmf = ACWT.makeqmfpair(wt)
	wp = zeros(n,L+1)
	wp[:,1] = x

    @inbounds begin
    	for d=0:(L-1)
    		for b=0:(2^d-1)
    		   s = wp[echant(n,d,b),1]
    		   high = acfilter(s,Qmf)
    		   low = acfilter(s,Pmf)
    		   wp[echant(n,d,b),L+1-d] = high
    		   wp[echant(n,d,b),1] = low
    		 end
    	end
    end
	return wp
end

"""Computes the column-wise forward acwat coeficients for 2D signals."""
function hacwt(x::AbstractArray{Float64,2},
                         wt::OrthoFilter,
                         L::Integer=maxtransformlevels(x[:,1]))
    nrow, ncol = size(x)
    W = Array{Float64,3}(undef,nrow,L+1,ncol)
    @inbounds begin
        for i in 1:ncol
            W[:,:,i] = acwt(x[:,i],wt,L)
        end
    end
    return W
end

"""Computes the row-wise 1D acwt coeficients for 2D signals."""
function vacwt(x::AbstractArray{Float64,2},
                      wt::OrthoFilter,
                      L::Integer=maxtransformlevels(x[1,:]))
    nrow, ncol = size(x)
    W = Array{Float64,3}(undef,ncol,L+1,nrow)
    @inbounds begin
        for i in 1:nrow
            W[:,:,i] = acwt(x[i,:],wt,L)
        end
    end
    return W
end

function acwt(x::AbstractArray{Float64,2}, wt::OrthoFilter,
              Lrow::Integer=maxtransformlevels(x[1,:]),
              Lcol::Integer=maxtransformlevels(x[:,1]))
    nrow, ncol = size(x)
    W3d = hacwt(x,wt,Lcol)
    W4d = Array{Float64,4}(undef,Lcol+1,ncol,Lrow+1,nrow)
    @inbounds begin
        for i in 1:Lcol+1
            W4d[i,:,:,:] = vacwt(W3d[:,i,:],wt,Lrow)
        end
    end
    W4d = permutedims(W4d, [4,2,3,1])
    return W4d
end

# Inverse Transforms
function iacwt(x::AbstractArray{Float64,4})
    nrow, ncol, Lrow, Lcol = size(x)
    W4d = permutedims(x,[4,2,3,1])
    W3d = Array{Float64,3}(undef, nrow, Lcol, ncol)
    @inbounds begin
        for i in 1:Lcol
            for j in 1:nrow
                W3d[j,i,:] = iacwt(W4d[i,:,:,j])
            end
        end
    end
    y = Array{Float64,2}(undef, nrow, ncol)
    @inbounds begin
        for i in 1:ncol
            y[:,i] = iacwt(W3d[:,:,i])
        end
    end
    return y
end

# Wavelet Packets functions
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

function acwpt(x::Vector{T},
               wt::OrthoFilter,
               L::Integer=maxtransformlevels(x)) where T<:Number

    root = AcwptNode(x)
    acwpt(x,root,wt,L)
    return root
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
                node.parent.right = AcwptNode(node.data,node)
            else
                node.parent.left = AcwptNode(node.data,node)
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

end # End module
