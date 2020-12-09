using
    AbstractTrees,
    AutocorrelationShell,
    BenchmarkTools,
    CSV,
    DataFrames,
    LinearAlgebra,
    Plots,
    Random,
    Wavelets

## Setup
Q = qfilter(wavelet(WT.db2));
P = pfilter(wavelet(WT.db2));

x = range(0,stop=1,length=2^8);
wavelet_test_data = CSV.File("./test/experiments/tidytree/wavelet_test_256.csv") |> DataFrame
rng = MersenneTwister(123);

## Functions
function collect_leaves(x::BinaryNode) return collect(Leaves(x)) end

function threshold_bestbasis!(x::BinaryNode, type::String, t::Float64)
    @inbounds begin
        for n in collect_leaves(x)
            acthreshold!(n.data, type, t)
        end
    end
end

function threshold_bestbasis(x::BinaryNode, type::String, t::Float64)
    y = deepcopy(x)
    threshold_bestbasis!(y, type, t)
    return y
end

function acwpt_snr(x::AbstractArray{T}, y::AbstractArray{T}, type::String, step::Float64) where T<:Number
    bb = acwpt(y, P, Q) |> acwptPostOrderBestBasis
    leaves = collect_leaves(bb)
    max_coef = collect(
        Iterators.flatten(map(x -> x.data, leaves))
    ) |> (x -> abs.(x)) |> maximum

    thresh = [0:step:round(max_coef);]
    reconst = map(t -> threshold_bestbasis(bb, type, t), thresh) |> (x -> map(iacwpt, x))
    _snr = [snr(x, r) for r in reconst]

    return thresh, _snr
end

function acwt_snr(y::AbstractArray, x::AbstractArray, type, step)
    decomp = acwt(x; L=0, P=P, Q=Q)
    max_coef = decomp |> (a -> abs.(a)) |> maximum
    thresh = [0:step:round(max_coef);]

    reconst = map(t -> acthreshold(decomp, type, t), thresh) |> (a -> map(iacwt, a))
    _snr = [snr(y, r) for r in reconst]

    return thresh, _snr
end

function test_denoise(y::AbstractArray, noise::Float64, type::String; rng=MersenneTwister(123))
    y_noisy = make_noisy(y, rng, noise);

    bb = acwpt(y_noisy, P, Q) |> acwptPostOrderBestBasis
    wt_decomp = acwt(y_noisy; L=1, P, Q);

    t, s = acwt_snr(y, y_noisy, type, 0.01);
    wt_thresh = acthreshold(wt_decomp, type, t[findmax(s)[2]]);
    wt_reconst = iacwt(wt_thresh);

    t, s = acwpt_snr(y, y_noisy, type, 0.01);
    wpt_thresh = threshold_bestbasis(bb, type, t[findmax(s)[2]]);
    wpt_reconst = iacwpt(wpt_thresh);

    base = norm(y-y_noisy)
    acwt_res = norm(y-wt_reconst)
    acwpt_res = norm(y-wpt_reconst)

    @printf(
        "Baseline: %.4f, ACWT: %.4f, ACWPT: %.4f\n",
        base,
        acwt_res,
        acwpt_res
    )

    return base, acwt_res, acwpt_res
end

## Quadchirp
y = wavelet_test_data.quadchirp;
y_noisy = make_noisy(y, rng, 0.07);
p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")

function tree_to_bitmap(x::BinaryNode)
    nrow, ncol = length(x.data), dyadlength(x.data) + 1
    arr, hash = zeros(nrow, ncol), zeros(ncol - 1)

    @inbounds begin
        for n in collect(PreOrderDFS(x))
            d = n.depth
            l = length(n.data)/2^d
            _start = Int(hash[d]+1)
            _end = Int(_start + l - 1)
            if !isdefined(n, :left)
                hash[d+1:end] .= _start + l/2 - 1
                arr[_start:_end, d+1] .= 1
            end
            if !isdefined(n, :right)
                hash[d+1:end] .= _end
                arr[_start:_end, d+1] .= 1
            end
            hash[d] += l
        end
    end

    return arr
end

function selectednodes_plot(x::BinaryNode; nodecolor::Symbol = :red)
    bitmap = tree_to_bitmap(x)
    nrow, ncol = size(bitmap)
    p = heatmap(
        transpose(bitmap),
        color = [:black, nodecolor],
        yflip=true,
        legend=false,
        xlims = (1, nrow+0.5),
        ylims = (0.5, ncol+0.5),
        xticks = false,
        yticks = 0:ncol
    )
    hline!([1.5:1:(ncol-0.5);], color=:white)
    @inbounds begin
        for i in 1:ncol
            for j in 1:2^(i-1)
                vpos = (nrow/2^i)*(2*j-1) + 0.5
                plot!(vpos*ones(ncol-i+1), (i+0.5):(ncol+0.5), color=:white)
            end
        end
    end
    return p
end

acwpt_decomp = acwpt(y_noisy, P, Q);
bb = acwptPostOrderBestBasis(acwpt_decomp);
bb_plot = selectednodes_plot(bb)

acwt_decomp = acwt(y_noisy; L=1, P, Q);

p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")
plot!(x, acwt_reconst, label = "ACWT Reconstructed")
plot!(x, acwpt_reconst, label = "ACWPT Reconstructed")

baseline = zeros(10);
acwt_score = zeros(10);
acwpt_score = zeros(10);

for i in 1:10
    baseline[i], acwt_score[i], acwpt_score[i] = test_denoise(
        wavelet_test_data.quadchirp,
        0.07,
        "hard";
        rng = MersenneTwister()
    )
end

mean(baseline)
mean(acwt_score)
mean(acwpt_score)

## Doppler
y = wavelet_test_data.doppler;
y_noisy = make_noisy(y, rng, 0.07);
p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")

acwpt_decomp = acwpt(y_noisy, P, Q);
bb = acwptPostOrderBestBasis(acwpt_decomp);
bb_plot = plot_tfbdry(bb)

acwt_decomp = acwt(y_noisy; L=1, P, Q);

t, s = acwt_snr(y, y_noisy, "soft", 0.01);
t = acthreshold(acwt_decomp, "soft", t[findmax(s)[2]]);
acwt_reconst = iacwt(t)

t, s = acwpt_snr(y, y_noisy, "soft", 0.01);
t = threshold_bestbasis(bb, "soft", t[findmax(s)[2]]);
acwpt_reconst = iacwpt(t);

p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")
plot!(x, acwt_reconst, label = "ACWT Reconstructed")
plot!(x, acwpt_reconst, label = "ACWPT Reconstructed")

norm(y - y_noisy)
norm(y - acwt_reconst)
norm(y - acwpt_reconst)

baseline = zeros(10);
acwt_score = zeros(10);
acwpt_score = zeros(10);

for i in 1:10
    baseline[i], acwt_score[i], acwpt_score[i] = test_denoise(
        wavelet_test_data.doppler,
        0.07,
        "hard";
        rng = MersenneTwister()
    )
end

mean(baseline)
mean(acwt_score)
mean(acwpt_score)

## Heavy Sine
y = wavelet_test_data.heavy_sine;
y_noisy = make_noisy(y, rng, 0.5);
p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")

acwpt_decomp = acwpt(y_noisy, P, Q)
bb = acwptPostOrderBestBasis(acwpt_decomp)
bb_plot = plot_tfbdry(bb)

acwt_decomp = acwt(y_noisy; L=1, P, Q);

t, s = acwt_snr(y, y_noisy, "soft", 0.01);
t = acthreshold(acwt_decomp, "soft", t[findmax(s)[2]]);
acwt_reconst = iacwt(t)

t, s = acwpt_snr(y, y_noisy, "soft", 0.01);
t = threshold_bestbasis(bb, "soft", t[findmax(s)[2]]);
acwpt_reconst = iacwpt(t);

p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")
plot!(x, acwt_reconst, label = "ACWT Reconstructed")
plot!(x, acwpt_reconst, label = "ACWPT Reconstructed")

baseline = zeros(10);
acwt_score = zeros(10);
acwpt_score = zeros(10);

for i in 1:10
    baseline[i], acwt_score[i], acwpt_score[i] = test_denoise(
        wavelet_test_data.heavy_sine,
        0.5,
        "hard";
        rng = MersenneTwister()
    )
end

mean(baseline)
mean(acwt_score)
mean(acwpt_score)

## Mish Mash
y = wavelet_test_data.mishmash;
y_noisy = make_noisy(y, rng, 0.5);
p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")

acwpt_decomp = acwpt(y_noisy, P, Q);
bb = acwptPostOrderBestBasis(acwpt_decomp);
bb_plot = plot_tfbdry(bb)

acwt_decomp = acwt(y_noisy; L=1, P, Q);

t, s = acwt_snr(y, y_noisy, "soft", 0.01);
t = acthreshold(acwt_decomp, "soft", t[findmax(s)[2]]);
acwt_reconst = iacwt(t)

t, s = acwpt_snr(y, y_noisy, "soft", 0.01);
t = threshold_bestbasis(bb, "soft", t[findmax(s)[2]]);
acwpt_reconst = iacwpt(t);

p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")
plot!(x, acwt_reconst, label = "ACWT Reconstructed")
plot!(x, acwpt_reconst, label = "ACWPT Reconstructed")

baseline = zeros(10);
acwt_score = zeros(10);
acwpt_score = zeros(10);

for i in 1:10
    baseline[i], acwt_score[i], acwpt_score[i] = test_denoise(
        wavelet_test_data.mishmash,
        0.5,
        "hard";
        rng = MersenneTwister()
    )
end

mean(baseline)
mean(acwt_score)
mean(acwpt_score)

## Block
y = wavelet_test_data.blocks;
y_noisy = make_noisy(y, rng, 0.5);
p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")

acwpt_decomp = acwpt(y_noisy, P, Q);
bb = acwptPostOrderBestBasis(acwpt_decomp);
bb_plot = plot_tfbdry(bb)

acwt_decomp = acwt(y_noisy; L=1, P, Q);

t, s = acwt_snr(y, y_noisy, "soft", 0.01);
t = acthreshold(acwt_decomp, "soft", t[findmax(s)[2]]);
acwt_reconst = iacwt(t)

t, s = acwpt_snr(y, y_noisy, "soft", 0.01);
t = threshold_bestbasis(bb, "soft", t[findmax(s)[2]]);
acwpt_reconst = iacwpt(t);

p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")
plot!(x, acwt_reconst, label = "ACWT Reconstructed")
plot!(x, acwpt_reconst, label = "ACWPT Reconstructed")

baseline = zeros(10);
acwt_score = zeros(10);
acwpt_score = zeros(10);

for i in 1:10
    baseline[i], acwt_score[i], acwpt_score[i] = test_denoise(
        wavelet_test_data.blocks,
        0.5,
        "hard";
        rng = MersenneTwister()
    )
end

mean(baseline)
mean(acwt_score)
mean(acwpt_score)

## Bumps
y = wavelet_test_data.bumps;
y_noisy = make_noisy(y, rng, 0.5);
p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")

acwpt_decomp = acwpt(y_noisy, P, Q);
bb = acwptPostOrderBestBasis(acwpt_decomp);
bb_plot = plot_tfbdry(bb)

acwt_decomp = acwt(y_noisy; L=1, P, Q);

t, s = acwt_snr(y, y_noisy, "soft", 0.01);
t = acthreshold(acwt_decomp, "soft", t[findmax(s)[2]]);
acwt_reconst = iacwt(t)

t, s = acwpt_snr(y, y_noisy, "soft", 0.01);
t = threshold_bestbasis(bb, "soft", t[findmax(s)[2]]);
acwpt_reconst = iacwpt(t);

p = plot(x, y, label = "Clean")
plot!(x, y_noisy, label = "Noisy")
plot!(x, acwt_reconst, label = "ACWT Reconstructed")
plot!(x, acwpt_reconst, label = "ACWPT Reconstructed")

baseline = zeros(10);
acwt_score = zeros(10);
acwpt_score = zeros(10);

for i in 1:10
    baseline[i], acwt_score[i], acwpt_score[i] = test_denoise(
        wavelet_test_data.bumps,
        0.5,
        "hard";
        rng = MersenneTwister()
    )
end

mean(baseline)
mean(acwt_score)
mean(acwpt_score)
