module ACUtil
export
    collectleaves,
    NormEntropy,
    findcost,
    makenoisy,
    snr,
    acwt_heatmap
using ..ACWT
using AbstractTrees, LinearAlgebra, Random, Wavelets, Plots

"""Collects the leaves of the best basis tree"""
collectleaves(x::AcwptNode) = collect(Leaves(x))

function collectleaves(coefvec::Vector{Integer}, bt::BitVector, i::Integer=1)
    M = length(bt)
    if rightchild(i) <= M
        if bt[leftchild(i)] != 0
            collectleaves(coefvec,bt,leftchild(i))
        end
        if bt[rightchild(i)] !=  0
            collectleaves(coefvec,bt,rightchild(i))
        end
        if bt[leftchild(i)] == 0 & bt[rightchild(i)] == 0
            append!(coefvec,i)
        end
    else
        if bt[i] == 1
            append!(coefvec,i)
        end
    end
end

# By its nature the indexes are in PostOrder traversal oreder
function collectleaves(bt::BitVector, exclude_scaling::Bool=false)
    coefvec = Vector{Integer}(undef,0)
    collectleaves(coefvec,bt,1)
    if exclude_scaling # excluding scaling coefficients
        popfirst!(coefvec)
    end
    return coefvec
end

## Thresholding
function Wavelets.Threshold.threshold!(x::AcwptNode, TH::Wavelets.Threshold.THType, t::Real)
    @assert t >= 0
    for leaf in collectleaves(x)
        Wavelets.Threshold.threshold!(leaf.data,TH,t)
    end
end

function Wavelets.Threshold.threshold(x::AcwptNode, TH::Wavelets.Threshold.THType, t::Real)
    y = deepcopy(x)
    Wavelets.Threshold.threshold!(y,TH,t)
    return y
end

## Entropy
struct NormEntropy <: Wavelets.Entropy end

function Wavelets.Threshold.coefentropy(x::T, et::NormEntropy, nrm::T, p::T=1) where T<:AbstractFloat
    s = abs((x/nrm))
    return s^p
end

function Wavelets.Threshold.coefentropy(x::AbstractArray{T}, et::NormEntropy,
                     nrm::T=norm(x),
                     p::Float64=1.0) where T<:AbstractFloat
    @assert nrm >= 0
    sum = zero(T)
    nrm == sum && return sum
    for i in eachindex(x)
        @inbounds sum += Wavelets.Threshold.coefentropy(x[i], et, nrm, p)
    end
    return sum
end

# Computes vector of cost values
"""
    findcost(W, et)

Find the cost of each node in the stationary wavelet matrix `W` using a given cost function `et`.
"""
function findcost(W::Array{Float64,2}, et::Wavelets.Entropy)
    return [Wavelets.Threshold.coefentropy(w,et) for w in eachcol(W)]
end

## Others
"""
    makenoisy(x::AbstractArray{<:Number}, rng::MersenneTwister, a::Float64)

Adds gaussian white noise to an image.

# Arguments
`x::AbstractArray{<:Number}`: image
`rng::MersenneTwister`: random number generator ex) rng = MersenneTwister(123)
`a::Float64`: A scaling parameter to adjust level of noise.

# Example
```julia
using AutocorrelationShell, Images, FileIO, Random

# Load test image
img = load("./test/pictures/boat.jpg")
img = Float64.(Gray.(img))

noisy_image = makenoisy(img, MersenneTwister(123), 0.7)
```
"""
function makenoisy(x::AbstractArray{<:Number}, rng::MersenneTwister, a::Float64)
    noisy = x + a * randn(rng, Float64, size(x))
    return noisy
end

"""
    snr(f::AbstractArray{T}, g::AbstractArray{T}) where T<:Number

Computes the signal to noise ratio between two signals. Returns signal to noise ratio(dB).
"""
function snr(f::AbstractArray{T}, g::AbstractArray{T}) where T<:Number
    L2_f = norm(f)
    L2_fg = norm(f-g)
    return 20*log10(L2_f/L2_fg)
end

"""
    acwt_heatmap(x::AbstractArray{<:Number,2})

Takes the natural log of the absolute value of each cell, and plots a heatmap.
Can be used for visualizing an input image or the coefficient matrix of a 2D ACW decomposition.

# Arguments
- `x::AbstractArray{<:Number,2}`: A image or a 2D matrix
"""
function acwt_heatmap(x::AbstractArray{<:Number,2})
    heatmap(log.(abs.(x)),
            yflip=true,
            axis=nothing,
            colorbar_entry=false,
            aspect_ratio=:equal,
            showaxis=false)
end

end # End module
