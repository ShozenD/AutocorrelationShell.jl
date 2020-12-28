module ACUtil
export
    collectleaves,
    threshold!,
    threshold,
    NormEntropy,
    coefentropy
using ..ACWT
using AbstractTrees, LinearAlgebra, Wavelets

"""Collects the leaves of the best basis tree"""
collectleaves(x::AcwptNode) = collect(Leaves(x))


## Thresholding
function Wavelets.Threshold.threshold!(x::AcwptNode, TH::Wavelets.Threshold.THType, t::Real)
    @assert t >= 0
    for leaf in collectleaves(x)
        threshold!(leaf.data,TH,t)
    end
end

function Wavelets.Threshold.threshold(x::AcwptNode, TH::Wavelets.Threshold.THType, t::Real)
    y = deepcopy(x)
    threshold!(y,TH,t)
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

## Others
"""
    snr(f::AbstractArray{T}, g::AbstractArray{T}) where T<:Number

Computes the signal to noise ratio between two signals. Returns signal to noise ratio(dB).
"""
function snr(f::AbstractArray{T}, g::AbstractArray{T}) where T<:Number
    L2_f = norm(f)
    L2_fg = norm(f-g)
    return 20*log10(L2_f/L2_fg)
end

end # End module
