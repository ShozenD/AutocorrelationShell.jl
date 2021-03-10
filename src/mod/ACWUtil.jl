module ACWUtil
using ..ACWT, ..ACTransforms
using AbstractTrees, Wavelets, Random

"""
    make_noisy(x::AbstractArray{<:Number}, rng::MersenneTwister, a::Float64)

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

rng = MersenneTwister(123)
noisy_image = make_noisy(img, rng, 0.7)
```
"""
function make_noisy(x::AbstractArray{<:Number}, rng::MersenneTwister, a::Float64)
    noisy = x + a * randn(rng, Float64, size(x))
    return noisy
end

"""
    get_snr(y::AbstractArray{T, Integer}, x::AbstractArray{T, Integer}; type::AbstractString="hard", step::Float64=0.5) where T<:Number

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
function get_snr(y::AbstractArray{T, Integer}, x::AbstractArray{T, Integer}; type::AbstractString="hard", step::Float64=0.5) where T<:Number
    max_coef = maximum(abs.(x)) # find largest coefficient
    num_coef = length(x)
    snr_list = zeros(0)
    coef_ratio_list = zeros(0)
    @inbounds begin
        for i in range(0, stop=round(max_coef), step=step)
            thresh = acthreshold(x, type, i)
            reconst = iac2d(thresh)
            num_nonzero_coef = sum(abs.(thresh) .> 0)
            append!(snr_list, snr(y, reconst))
            append!(coef_ratio_list, num_nonzero_coef/num_coef)
            println(i/round(max_coef) * 100)
        end
    end
    return coef_ratio_list, snr_list
end

"""
    get_ssim(y::AbstractArray{T}, x::AbstractArray{T}; type::AbstractString="hard", step::Float64=0.5) where T<:Number

Thresholds the wavelet coefficients matrix of a noisy signal, reconstructs it, then
computes the structural similarity index between the reconstructed signal and the
original non-noisy signal. Returns a list of the ratio of non-zero coefficient at
each threshold and a list containing the SSIM values.

# Arguments
- `y`: original 2D signal
- `x`: wavelet coefficient matrix of a noisy 2D signal
- `type`: type of thresholding (soft or hard)
- `step`: step at which to increase the threshold

# Example
```julia
using AutocorrelationShell, Wavelets, FileIO, Images, ImageQualityIndexes, Random

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

ac_noisy = acwt2d(noisy; L=2, P=P, Q=Q)

coef_ratios, ssim_values = get_ssim(img, ac_noisy; type="hard", step=0.5)
```
"""
function get_ssim(y::AbstractArray{T}, x::AbstractArray{T}; type::AbstractString="hard", step::Float64=0.5) where T<:Number
    max_coef = maximum(abs.(x)) # find largest coefficient
    num_coef = length(x)
    ssim_list = zeros(0)
    coef_ratio_list = zeros(0)
    @inbounds begin
        for i in range(0, stop=round(max_coef), step=step)
            thresh = acthreshold(x, type, i)
            reconst = iac2d(thresh)
            num_nonzero_coef = sum(abs.(thresh) .> 0)
            append!(ssim_list, assess(SSIM(), y, reconst))
            append!(coef_ratio_list, num_nonzero_coef/num_coef)
            println(i/round(max_coef) * 100)
        end
    end
    return coef_ratio_list, ssim_list
end

"""
    acwt_heatmap(x::AbstractArray{<:Number})

Takes the natural log of the absolute value of each cell, and plots a heatmap.
Can be used for visualizing an input image or the coefficient matrix of a 2D ACW decomposition.

# Arguments
- `x::AbstractArray{<:Number}`: A image or a 2D matrix
"""
function acwt_heatmap(x::AbstractArray{<:Number})
    heatmap(log.(abs.(x)),
            yflip=true,
            axis=nothing,
            colorbar_entry=false,
            aspect_ratio=:equal,
            showaxis=false)
end

"""
    acwt_snr(y::AbstractArray, x::AbstractArray, type, step)

Calculates the signal to noise ratio of a given signal and a noisy version of the same signal.
"""
function acwt_snr(y::AbstractArray, x::AbstractArray, type, step)
    decomp = acwt(x; L=0, P=pfilter(wavelet(WT.db2)), Q=qfilter(wavelet(WT.db2)))
    max_coef = decomp |> (a -> abs.(a)) |> maximum
    thresh = [0:step:round(max_coef);]

    reconst = map(t -> acthreshold(decomp, type, t), thresh) |> (a -> map(iacwt, a))
    _snr = [snr(y, r) for r in reconst]

    return thresh, _snr
end

"""
    acwpt_snr(x::AbstractArray{T}, y::AbstractArray{T}, type::String, step::Float64) where T<:Number

Returns the threshold values and corresponding SNR values
"""
function acwpt_snr(x::AbstractArray{T}, y::AbstractArray{T}, type::String, step::Float64) where T<:Number
    bb = acwpt(y, pfilter(wavelet(WT.db2)), qfilter(wavelet(WT.db2))) |> acwpt__bb
    leaves = collect_leaves(bb)
    max_coef = collect(
        Iterators.flatten(map(x -> x.data, leaves))
    ) |> (x -> abs.(x)) |> maximum

    thresh = [0:step:round(max_coef);]
    reconst = map(t -> threshold_bestbasis(bb, type, t), thresh) |> (x -> map(iacwpt, x))
    _snr = [snr(x, r) for r in reconst]

    return thresh, _snr
end


end
