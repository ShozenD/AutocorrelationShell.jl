"""
    make_noisy(x, rng, a)

Adds gaussian white noise to an image.

# Arguments
`x::AbstractArray{<:Number}`: image
`rng`: random number generator ex) rng = MersenneTwister(123)
`a::Float64`: A scaling parameter to adjust level of noise.
"""
function make_noisy(x::AbstractArray{<:Number}, rng, a::Float64)
    noisy = x + a * randn(rng, Float64, size(x))
    return noisy
end

"""
    snr(f::AbstractArray{T}, g::AbstractArray{T}) where T<:Number

Computes the signal to noise ratio between two signals. Returns signal to noise ratio(dB).

# Arguments
`f::AbstractArray{<:Number}`: Original signal
`g::AbstractArray{<:Number}`: Noisy signal
"""
function snr(f::AbstractArray{T}, g::AbstractArray{T}) where T<:Number
    L2_f = norm(f)
    L2_fg = norm(f-g)
    return 20*log10(L2_f/L2_fg)
end

"""
    get_snr(y, x; type, step)

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

coef_ratio, snr_list = get_snr(img, ac_noisy, type="soft", step=0.5)
```
"""
function get_snr(y::AbstractArray{T}, x::AbstractArray{T}; type::AbstractString="hard", step::Float64=0.5) where T<:Number
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

## Entropy functions
# ENTROPY TYPES

# Extend Wavetes.Entropy to include norm and threshold entropy methods
struct NormEntropy <: Wavelets.Entropy end
struct ThresholdEntropy <: Wavelets.Entropy end # Should we implement?

# Entropy measures: Additive with wentropy(0) = 0
# all coefs assumed to be on [-1,1] after normalization with nrm
# given x and y, where x has "more concentrated energy" than y
# then coefentropy(x, et, norm) <= coefentropy(y, et, norm) should be satisfied.

# (non-normalized) Shannon Entropy
function wentropy(x::T, et::Wavelets.ShannonEntropy, nrm::T)  where T<:AbstractFloat
    s = (x/nrm)^2
    if s == 0.0
        return -zero(T)
    else
        return -s*log(s)
    end
end

# log energy entropy
function wentropy(x::T, et::Wavelets.LogEnergyEntropy, nrm::T) where T<:AbstractFloat
    s = (x/nrm)^2
    if s == 0.0
        return -zero(T)
    else
        return -log(s)
    end
end

function wentropy(x::T, et::NormEntropy, nrm::T; p::T=1) where T<:AbstractFloat
    s = (x/nrm)^2
    return s^p
end

function wentropy(x::AbstractArray{T}, et::Wavelets.Entropy, nrm::T=norm(x)) where T<:AbstractFloat
    @assert nrm >= 0
    sum = zero(T)
    nrm == sum && return sum
    for i in eachindex(x)
        @inbounds sum += wentropy(x[i], et, nrm)
    end
    return sum
end
