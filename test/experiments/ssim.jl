using Plots, Wavelets, LinearAlgebra, FileIO, Images, Random
include("../src/AutocorrelationShell.jl")
using Main.AutocorrelationShell

# Load test image
img = load("./test/pictures/lenna.jpg")
img = Float64.(Gray.(img))

rng = MersenneTwister(123)
noisy = make_noisy(img, rng, 0.7)

# Access image quality using SSIM
using ImageQualityIndexes

assess(SSIM(), img, noisy)

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

H = wavelet(WT.db2)
L = 4
Q = qfilter(H)
P = pfilter(H)

ac_noisy = ac2d(noisy, L, P, Q)

function get_ssim(y, x, type, step)
    max_coef = maximum(abs.(x)) # find largest coefficient
    num_coef = length(x)
    ssim_list = zeros(0)
    coef_ratio_list = zeros(0)
    for i in range(0, stop=round(max_coef), step=step)
        thresh = acthreshold(x, type, i)
        reconst = iac2d(thresh)
        num_nonzero_coef = sum(abs.(thresh) .> 0)
        append!(ssim_list, assess(SSIM(), y, reconst))
        append!(coef_ratio_list, num_nonzero_coef/num_coef)
        println(i/round(max_coef) * 100)
    end
    return coef_ratio_list, ssim_list
end

coef_ratio_list, ssim_list = get_ssim(img, ac_noisy, "soft", 0.2)

max_coef = maximum(abs.(ac_noisy))
p1 = plot(range(0, stop=round(max_coef), step=0.2), ssim_list, xlabel="Threshold", ylabel="SSIM")
p2 = plot(coef_ratio_list, ssim_list, xlabel="Ratio of Non-zero Coefficients", ylabel="SSIM")
p = plot(p1, p2, layout = (2), legend=false)
