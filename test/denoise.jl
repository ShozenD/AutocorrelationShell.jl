# --- Test denoising algorithm ---
# Import packages
using Plots, Wavelets, LinearAlgebra, FileIO, Images, Random

include("../src/AutocorrelationShell.jl")
using Main.AutocorrelationShell

# Load test image
img = load("./test/pictures/lenna.jpg")
img = Float64.(Gray.(img))

rng = MersenneTwister(123)
noisy = make_noisy(img, rng, 0.7)

# Apply Wavelet transform
H = wavelet(WT.db2)
L = 4
Q = qfilter(H)
P = pfilter(H)

ac_noisy = ac2d(noisy, L, L, P, Q)
acwt2D(noisy; P=P, Q=Q)
iacwt2D(ac_noisy)

heatmap(ac_noisy[:,:,5,5])

# Thresholding
coef_ratio, snr_list = get_snr(img, ac_noisy, "soft", 0.2)

# Visualization
max_coef = maximum(abs.(ac_noisy))
p1 = plot(range(0, stop=round(max_coef), step=0.2), snr_list, xlabel="Threshold", ylabel="SNR(Db)")
p2 = plot(coef_ratio, snr_list, xlabel="Ratio of Non-zero Coefficients", ylabel="SNR(Db)")
p = plot(p1, p2, layout = (2), legend=false)

savefig(p, "SNR")

findmax(snr_list)
range(0, stop=round(max_coef), step=0.2)[7]
coef_ratio[7]

denoised = acthreshold(ac_noisy, "soft", 1.20)
reconst = iac2d(denoised)

# Create Heatmap
h1 = heatmap(img, title="Original", legend=false)
h2 = heatmap(noisy, title="Noisy", legend=false)
h3 = heatmap(reconst, title="Denoised", legend=false)
h = plot(h1, h2, h3, layout=(1, 3), legend=false, size=(1200, 400))

savefig(h, "OND")
