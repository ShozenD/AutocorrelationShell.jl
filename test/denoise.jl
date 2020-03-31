# --- Test denoising algorithm ---
# Import packages
using Plots, Wavelets, LinearAlgebra, FileIO, Images
using Random

include("../src/AutocorrelationShell.jl")
using Main.AutocorrelationShell

# Load test image
img = load("./test/pictures/lenna.jpg")
img = Float64.(Gray.(img))

rng = MersenneTwister(123)

noisy = make_noisy(img, rng, 0.7)

# Apply Wavelet transform
H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

ac_noisy = ac2d(noisy, L, P, Q)

max_coef = maximum(abs.(ac2d_flatten(ac_noisy)))

# --- Thresholding ---
function get_snr(y, x, type, step)
    max_coef = maximum(abs.(ac2d_flatten(x))) # find largest coefficient
    snr_list = zeros(0)
    for i in range(0, stop=round(max_coef), step=step)
        thresh = acthreshold(x, type, i)
        reconst = iac2d(thresh)
        append!(snr_list, snr(y, reconst))
        println(i/round(max_coef) * 100)
    end
    return snr_list
end

function get_snr_ratio(y, x, type, step)
    max_coef = maximum(abs.(ac2d_flatten(x)))
    snr_list = zeros(0)
    coef_ratio_list = zeros(0)
    num_coef = length(ac2d_flatten(x))

    print(num_coef)

    for i in range(0, stop=round(max_coef), step=step)
        thresh = acthreshold(x, type, i)
        reconst = iac2d(thresh)
        num_nonzero_coef = sum(abs.(ac2d_flatten(thresh)) .> 0)
        append!(coef_ratio_list, num_nonzero_coef/num_coef)
        append!(snr_list, snr(y, reconst))
        println(i/round(max_coef) * 100)
    end
    return coef_ratio_list, snr_list
end

result = get_snr(img, ac_noisy, "soft", 0.2)

ratio, snr_val = get_snr_ratio(img, ac_noisy, "soft", 0.2)

plot(range(0, stop=round(max_coef), step=0.2), result, label="SNR", lw=1.5, xlims=(0, 10))

plot(ratio, snr_val, label="SNR", lw=1.5)
plot!(title = "SNR vs Threshold", xlabel = "Threshold", ylabel = "SNR(db)")
savefig("SNR_lenna_soft_step=0.2.png")

plot(range(0, stop=round(max_coef), step=0.2), ratio)

inv_theta = zeros(0)
for i in range(0, stop=round(max_coef), step=0.2)
    append!(inv_theta, 1/i)
end
plot(inv_theta, result)
plot(ratio, result)

using BenchmarkTools

@benchmark iac2d(ac_noisy)


Array{Float64,3}(undef, 2, 2, 2)
