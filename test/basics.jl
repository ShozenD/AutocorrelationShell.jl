using Plots, Wavelets, FileIO, Images, LinearAlgebra

include("../src/AutocorrelationShell.jl")
using Main.AutocorrelationShell

H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

# simulate signal
x = zeros(256)
x[128] = 1

decomp = fwt_ac(x,L,P,Q)
wiggle(decomp, Overlap = false)

# -- 2D Autocorrelation Wavelet Transform
img = load("./test/pictures/lenna.jpg")
img = Float64.(Gray.(img))

decomp = ac2d(img, L, P, Q)

h1 = heatmap(img, title="Original" ,legend=false)
h2 = heatmap(decomp[:,:,5,5], title="Transform" ,legend=false)
h = plot(h1, h2, size=(800, 400))

savefig(h, "2d-decomp")
