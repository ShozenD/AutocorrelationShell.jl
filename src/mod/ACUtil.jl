module ACUtil
export
    # modify ac2d output
    ac2d_flatten,
    # modify image
    make_noisy,
    # signal to noise ratio
    snr
using Random
using LinearAlgebra

"""
    make_noisy(x, rng, a)

Adds gaussian white noise to an image.
### Arguments
`x`: image
`rng`: random number generator ex) rng = MersenneTwister(123)
`a`: A scaling parameter to adjust level of noise.
"""
function make_noisy(x, rng, a::Real)
    noisy = x + a * randn(rng, Float64, size(x))
    return noisy
end

"""
    snr(f::AbstractArray{<:Number}, g::AbstractArray{<:Number})

Computes the signal to noise ratio between two signals.
### Arguments
`f`: Original signal
`g`: Noisy signal
### Returns
signal to noise ratio(dB)
"""
function snr(f::AbstractArray{<:Number}, g::AbstractArray{<:Number})
    L2_f = norm(f)
    L2_fg = norm(f-g)
    return 20*log10(L2_f/L2_fg)
end

end # module
