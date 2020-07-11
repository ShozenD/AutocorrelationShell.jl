using Plots, Wavelets, LinearAlgebra, Statistics, Random, FileIO
using Images, ImageQualityIndexes
include("../src/AutocorrelationShell.jl")
using Main.AutocorrelationShell

"""
    plot_heatmap(x)

Plots a heatmap of an image.
### Arguments
`x`: image
"""
function plot_heatmap(x)
    heatmap(log.(abs.(x)),
            yflip=true,
            axis=nothing,
            colorbar_entry=false,
            aspect_ratio=:equal,
            showaxis=false)
end

"""
    threshold_ssim(img,L,P,Q,noise,step_size,type)

Thresholds coefficients of a noisy image using soft
or hard thresholding at various threshold levels.
### Arguments
`img`: image
`L`: decomposition level
`P`: low AC shell filter
`Q`: high AC shell filter
`noise`: level of noise
`step_size`: step at which to increase the threshold
`type`: type of thresholding ("soft" or "hard")
### Returns
An 2-element array where the first element is an array of
thresholds and the second element is an array of signal
to noise ratios.
"""
function threshold_psnr(img,L,P,Q,noise,step_size,type)
    rng = MersenneTwister(123)
    img_noise = make_noisy(img,rng,noise)
    ac = ac2d(img_noise,L,P,Q)

    ac_max = maximum(abs.(ac[:]))
    step = collect(0:step_size:ac_max)

    psnr_list = zeros(0)
    for s in step
        ac_thresh = acthreshold(ac, type, s)
        thresh_img = iac2d(ac_thresh)
        ac_psnr = assess(PSNR(), img, thresh_img)
        append!(psnr_list, ac_psnr)
    end

    return step, psnr_list
end

"""
    threshold_ssim_plot(img,low_noise,high_noise,step)

Plots the threshold against the signal to noise ratio for
    - an image with low noise using soft thresholding
    - an image with low noise using hard thresholding
    - an image with high noise using soft thresholding
    - an image with high noise using hard thresholding
This plot does not include labels or a legend.
### Arguments
`img`: image
`low_noise`: level of low noise
`high_noise`: level of high noise
`step`: step at which to increase the threshold
"""
function threshold_psnr_plot(img,low_noise,high_noise,step)
    x_label = "Threshold"
    y_label = "PSNR"

    a, b = threshold_psnr(img, L, P, Q, low_noise, step, "soft")
    c, d = threshold_psnr(img, L, P, Q, low_noise, step, "hard")
    e, f = threshold_psnr(img, L, P, Q, high_noise, step, "soft")
    g, h = threshold_psnr(img, L, P, Q, high_noise, step, "hard")

    plot(a, b, xlabel=x_label, ylabel=y_label, label="", legend=false,
        linewidth=2, color=:blue, linestyle=:dot)
    plot!(c, d, xlabel=x_label, ylabel=y_label,
        label="",legend=false,
        linewidth=2,color=:blue)
    plot!(e, f, xlabel=x_label, ylabel=y_label,
        label="",legend=false,
        linewidth=2,color=:red,linestyle=:dot)
    p = plot!(g, h, xlabel=x_label, ylabel=y_label,
        label="",legend=false,
        linewidth=2,color=:red)
    return p
end

"""
    threshold_snr_plot_leg(img,low_noise,high_noise,step)

Plots the threshold against the signal to noise ratio for
    - an image with low noise using soft thresholding
    - an image with low noise using hard thresholding
    - an image with high noise using soft thresholding
    - an image with high noise using hard thresholding
This plot includes a legend.
### Arguments
`img`: image
`low_noise`: level of low noise
`high_noise`: level of high noise
`step`: step at which to increase the threshold
"""
function threshold_psnr_plot_leg(img,low_noise,high_noise,step)
    x_label = "Threshold"
    y_label = "PSNR"

    a, b = threshold_psnr(img,L,P,Q,low_noise,step,"soft")
    c, d = threshold_psnr(img,L,P,Q,low_noise,step,"hard")
    e, f = threshold_psnr(img,L,P,Q,high_noise,step,"soft")
    g, h = threshold_psnr(img,L,P,Q,high_noise,step,"hard")
    plot(a,b,xlabel=x_label,ylabel=y_label,
        label="soft, low noise",legend=false,
        linewidth=2,color=:blue,linestyle=:dot)
    plot!(c,d,xlabel=x_label,ylabel=y_label,
        label="hard, low noise",legend=false,
        linewidth=2,color=:blue)
    plot!(e,f,xlabel=x_label,ylabel=y_label,
        label="soft, high noise",legend=false,
        linewidth=2,color=:red,linestyle=:dot)
    p=plot!(g,h,xlabel=x_label,ylabel=y_label,
        label="hard, high noise",legend=true,
        linewidth=2,color=:red)
    return p
end

# Generate threshold vs. snr plot for all images
PATH = "./test/pictures"
fn_list = [fn for fn in readdir(PATH) if fn != ".DS_Store"]
load_list = [load(string(PATH, "/", fn_list[i])) for i=1:length(fn_list)]
img_list = [Float32.(Gray.(load_list[i])) for i=1:length(load_list)]

H = wavelet(WT.db2)
L = 6
Q = qfilter(H)
P = pfilter(H)
plot_list = [threshold_psnr_plot(img_list[i],0.1,0.7,0.1) for i=1:length(img_list)]
airplane,baboon,barbara,boat,fruits,girlface,goldhill,lenna,peppers,wave = plot_list

# Create a figure of subplots of threshold vs. snr plots
using Plots.PlotMeasures
font1 = Plots.font("Helvetica", 12)
font2 = Plots.font("Helvetica", 7)
airplane_leg = threshold_psnr_plot_leg(img_list[1], 0.1, 0.7, 0.1)
plot_titles=["Airplane" "Baboon" "Barbara" "Boat" "Fruits" "Girlface" "Goldhill" "Lenna" "Peppers" "Wave"]

plotly()
plot(airplane_leg, baboon, barbara, boat, fruits,
    girlface,goldhill,lenna,peppers,wave,
    layout=(2,5),
    title=plot_titles,
    size=(1200,500),
    legend=true,
    titlefont=font1,
    guidefont=font2,
    xtickfont=font2,
    ytickfont=font2,
    bottom_margin = 10mm)
gui()

# Create a figure of subplots of the original images
plot(plot_heatmap(img_list[1]),
    plot_heatmap(img_list[2]),
    plot_heatmap(img_list[3]),
    plot_heatmap(img_list[4]),
    plot_heatmap(img_list[5]),
    plot_heatmap(img_list[6]),
    plot_heatmap(img_list[7]),
    plot_heatmap(img_list[8]),
    plot_heatmap(img_list[9]),
    plot_heatmap(img_list[10]),
    size=(1200,500),
    layout=(2,5))
gui()
