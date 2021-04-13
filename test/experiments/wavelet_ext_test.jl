using LinearAlgebra, Wavelets, AutocorrelationShell, WaveletsExt, Plots, DataFrames, CSV

PATH = "/Users/shozendan/Documents/autocorrelation-shell/test/data/wavelet_test_256.csv";
testdata = CSV.read(PATH, DataFrame);

x = testdata.doppler
x = x + 0.05 * randn(length(x))

y = swpd(x, wavelet(WT.db4))

y = acwpt(x, wavelet(WT.db4))

z = denoise(y, :swpd, nothing)

r = iswpd(z)

plot(x, label = "noisy")
plot!(r, label = "denoised")

