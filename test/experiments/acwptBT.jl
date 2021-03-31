using
    AutocorrelationShell,
    CSV,
    DataFrames,
    Random,
    LinearAlgebra,
    Statistics,
    Plots,
    Wavelets,
    ProgressBars,
    BenchmarkTools

include("./denoiseBestBasis.jl");

## Load Data
PATH = "/Users/shozendan/Documents/autocorrelation-shell/test/data/wavelet_test_256.csv";
testdata = CSV.read(PATH, DataFrame);

function testdenoise(x::Vector{Float64}, noise::Float64;
                     TH::Wavelets.Threshold.THType=HardTH(), N::Integer=50, M::Integer=100,
                     method::Symbol=:shozen)
    L₁ = Array{Float64,2}(undef,(N,3))
    L₂ = Array{Float64,2}(undef,(N,3))
    SNR = Array{Float64,2}(undef,(N,3))
    PSNR = Array{Float64,2}(undef,(N,3))

    X₀ = shiftsignal(x,M);
    wt = wavelet(WT.db4);
    @inbounds begin
        for n in ProgressBar(1:N)
            X₁, Y₀, Y₁, th = Σ(X₀, wt, n, noise, TH=TH, method=method)
            L₁[n,1], L₂[n,1], SNR[n,1], PSNR[n,1] = getresults(X₀,Y₀,Y₁,TH,th)
            L₁[n,2], L₂[n,2], SNR[n,2], PSNR[n,2] = getresults(X₀,Y₀,Y₁,TH,mean(th))
            L₁[n,3], L₂[n,3], SNR[n,3], PSNR[n,3] = getresults(X₀,Y₀,Y₁,TH,median(th))
        end
    end

    return mean(L₁,dims=1), mean(L₂,dims=1), mean(SNR,dims=1), mean(PSNR,dims=1)
end

function testdenoise_noshift(x::Vector{Float64}, noise::Float64;
                     TH::Wavelets.Threshold.THType=HardTH(), N::Integer=50, M::Integer=100,
                     method::Symbol=:shozen)
    L₁ = Array{Float64,2}(undef,(N,3))
    L₂ = Array{Float64,2}(undef,(N,3))
    SNR = Array{Float64,2}(undef,(N,3))
    PSNR = Array{Float64,2}(undef,(N,3))

    X₀ = noshiftsignal(x,M);
    wt = wavelet(WT.db4);
    @inbounds begin
        for n in ProgressBar(1:N)
            X₁, Y₀, Y₁, th = Σ(X₀, wt, n, noise, TH=TH, method=method)
            L₁[n,1], L₂[n,1], SNR[n,1], PSNR[n,1] = getresults(X₀,Y₀,Y₁,TH,th)
            L₁[n,2], L₂[n,2], SNR[n,2], PSNR[n,2] = getresults(X₀,Y₀,Y₁,TH,mean(th))
            L₁[n,3], L₂[n,3], SNR[n,3], PSNR[n,3] = getresults(X₀,Y₀,Y₁,TH,median(th))
        end
    end

    return mean(L₁,dims=1), mean(L₂,dims=1), mean(SNR,dims=1), mean(PSNR,dims=1)
end

## Shifted experiment
L1, L2, SNR, PSNR = testdenoise(testdata.bumps, 0.5, TH=SoftTH(), method=:jeff);

## Non-shifted experiment
L1, L2, SNR, PSNR = testdenoise_noshift(testdata.mishmash, 0.7, TH=SoftTH(), method=:jeff);

## Export results
df = DataFrame(L1 = vec(L1), L2 = vec(L2), PSNR = vec(PSNR), SNR = vec(SNR))
CSV.write("../../Desktop/results.csv", df)