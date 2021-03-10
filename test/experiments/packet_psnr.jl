using
    AutocorrelationShell,
    AbstractTrees,
    Wavelets,
    CSV,
    DataFrames,
    Random,
    Plotly,
    LinearAlgebra,
    Images,
    ImageQualityIndexes

## Setup
rng = MersenneTwister(123);

Q = qfilter(wavelet(WT.db2));
P = pfilter(wavelet(WT.db2));

wavelet_test_data = CSV.File(
        "/Users/shozendan/Documents/autocorrelation-shell/test/data/wavelet_test_256.csv"
    ) |> DataFrame;

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

function acwt_psnr(y::AbstractArray, x::AbstractArray, type, step)
    decomp = acwt(x; L=0, P=pfilter(wavelet(WT.db2)), Q=qfilter(wavelet(WT.db2)))
    max_coef = decomp |> (a -> abs.(a)) |> maximum
    thresh = [0:step:round(max_coef);]

    reconst = map(t -> acthreshold(decomp, type, t), thresh) |> (a -> map(iacwt, a))
    psnr = [assess(PSNR(), y, r) for r in reconst]

    return thresh, psnr
end

function acwpt_psnr(x::AbstractArray{T}, y::AbstractArray{T}, type::String, step::Float64) where T<:Number
    bb = acwpt(y, P, Q) |> acwpt_postorder_bb
    leaves = collect_leaves(bb)
    max_coef = collect(
        Iterators.flatten(map(x -> x.data, leaves))
    ) |> (x -> abs.(x)) |> maximum

    thresh = [0:step:round(max_coef);]
    reconst = map(t -> threshold_bestbasis(bb, type, t), thresh) |> (x -> map(iacwpt, x))
    psnr = [assess(PSNR(), x, r) for r in reconst]

    return thresh, psnr
end

y = wavelet_test_data.doppler;
decomp = acwpt(y, P, Q)
bb = acwpt_postorder_bb(decomp)
selectednodes_plot(bb)

## Doppler
y = wavelet_test_data.doppler;
y_noisy = make_noisy(y, rng, 0.07);

acw_th, acw_psnr = acwt_psnr(y, y_noisy, "soft", 0.01);
th, psnr = acwpt_psnr(y, y_noisy, "soft", 0.01);
Plotly.plot(acw_th, acw_psnr)
Plotly.plot(th, psnr)

maximum(acw_psnr)
maximum(psnr)

## Quadchirp
y = wavelet_test_data.quadchirp;
y_noisy = make_noisy(y, rng, 0.07);

acw_th, acw_psnr = acwt_psnr(y, y_noisy, "soft", 0.01);
th, psnr = acwpt_psnr(y, y_noisy, "soft", 0.01);
Plotly.plot(acw_th, acw_psnr)
Plotly.plot(th, psnr)

maximum(acw_psnr)
maximum(psnr)

## Heavy Sine
y = wavelet_test_data.heavy_sine;
y_noisy = make_noisy(y, rng, 0.5);

acw_th, acw_psnr = acwt_psnr(y, y_noisy, "soft", 0.01);
th, psnr = acwpt_psnr(y, y_noisy, "soft", 0.01);
Plotly.plot(acw_th, acw_psnr)
Plotly.plot!(th, psnr)

maximum(acw_psnr)
maximum(psnr)

## Mish Mash
y = wavelet_test_data.mishmash;
y_noisy = make_noisy(y, rng, 0.5);

acw_th, acw_psnr = acwt_psnr(y, y_noisy, "soft", 0.01);
th, psnr = acwpt_psnr(y, y_noisy, "soft", 0.01);
Plotly.plot(acw_th, acw_psnr)
Plotly.plot(th, psnr)

maximum(acw_psnr)
maximum(psnr)

## Block
y = wavelet_test_data.mishmash;
y_noisy = make_noisy(y, rng, 0.5);

acw_th, acw_psnr = acwt_psnr(y, y_noisy, "soft", 0.01);
th, psnr = acwpt_psnr(y, y_noisy, "soft", 0.01);
Plotly.plot(acw_th, acw_psnr)
Plotly.plot(th, psnr)

maximum(acw_psnr)
maximum(psnr)

## Bumps
y = wavelet_test_data.bumps;
y_noisy = make_noisy(y, rng, 0.5);

acw_th, acw_psnr = acwt_psnr(y, y_noisy, "soft", 0.01);
th, psnr = acwpt_psnr(y, y_noisy, "soft", 0.01);
Plotly.plot(acw_th, acw_psnr)
Plotly.plot(th, psnr)

maximum(acw_psnr)
maximum(psnr)
