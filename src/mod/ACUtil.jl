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

"""
    get_snr(y, x, type, step)

Thresholds the wavelet coefficients matrix of a noisy signal, reconstructs it, then
computes the signal to noise ratio between the reconstructed signal and the
original non-noisy signal. Returns a list of the ratio of non-zero coefficient at
each threshold and a list containing the signal to noise ratio.

# Arguments
- `y`: original 2D signal
- `x`: wavelet coefficient matrix of a noisy 2D signal
- `type`: type of thresholding (soft or hard)
- `step`: step at which to increase the threshold

# Example
```julia
img = load("./test/pictures/lenna.jpg")
img = Float64.(Gray.(img))

# Add noise to image
rng = MersenneTwister(123)
noisy = make_noisy(img, rng, 0.7)

# Apply Wavelet transform
H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

ac_noisy = ac2d(noisy, L, P, Q)

coef_ratio, snr_list = get_snr(img, ac_noisy, "soft", 0.5)
```
"""
function get_snr(y, x, type, step)
    max_coef = maximum(abs.(x)) # find largest coefficient
    num_coef = length(x)
    snr_list = zeros(0)
    coef_ratio_list = zeros(0)
    for i in range(0, stop=round(max_coef), step=step)
        thresh = acthreshold(x, type, i)
        reconst = iac2d(thresh)
        num_nonzero_coef = sum(abs.(thresh) .> 0)
        append!(snr_list, snr(y, reconst))
        append!(coef_ratio_list, num_nonzero_coef/num_coef)
        println(i/round(max_coef) * 100)
    end
    return coef_ratio_list, snr_list
end

end # module
