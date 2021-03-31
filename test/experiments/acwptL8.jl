using
    AutocorrelationShell,
    AbstractTrees,
    CSV,
    DataFrames,
    Random,
    LinearAlgebra,
    Plots,
    Statistics,
    Wavelets,
    ProgressBars

include("./denoiseBestBasis.jl")

## Load data
PATH = "/Users/shozendan/Documents/autocorrelation-shell/test/data/wavelet_test_256.csv";
testdata = CSV.read(PATH, DataFrame)

function maketrees(n, m, L) # n: length of original signal, m: number of signals
    y = BitArray{2}(undef,(n<<1-1,m))
    for idx in 1:m
        @inbounds y[:,idx] = maketree(n<<1,L+1)
    end
    return y
end

function L8Σ(X::Array{Float64,2}, wt::OrthoFilter, n::Integer, noise::Float64;
              TH::Wavelets.Threshold.THType=HardTH(), L::Integer=maxtransformlevels(X[:,1]),
              method::Symbol=:jeff)
    # add noise
    X₁ = addnoise(X,n,noise)
    # acwpt decompositions
    Y₀ = decompositions(X₁,wt,L)
    # Find best basis trees for each signal
    y₁ = maketrees(size(Y₀,1),size(Y₀,3),8)
    # Find threshold
    if TH == SoftTH()
        th = bestthreshold(Y₀,y₁,TH=TH,method=method,elbows=3)
    else
        th = bestthreshold(Y₀,y₁,TH=TH,method=method,elbows=2)
    end

    return X₁, Y₀, y₁, th
end

function L8testdenoise(x::Vector{Float64}, noise::Float64;
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
            X₁, Y₀, Y₁, th = L8Σ(X₀, wt, n, noise, TH=TH, method=method)
            L₁[n,1], L₂[n,1], SNR[n,1], PSNR[n,1] = getresults(X₀,Y₀,Y₁,TH,th)
            L₁[n,2], L₂[n,2], SNR[n,2], PSNR[n,2] = getresults(X₀,Y₀,Y₁,TH,mean(th))
            L₁[n,3], L₂[n,3], SNR[n,3], PSNR[n,3] = getresults(X₀,Y₀,Y₁,TH,median(th))
        end
    end

    return mean(L₁,dims=1), mean(L₂,dims=1), mean(SNR,dims=1), mean(PSNR,dims=1)
end

function L8testdenoise_noshift(x::Vector{Float64}, noise::Float64;
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
            X₁, Y₀, Y₁, th = L8Σ(X₀, wt, n, noise, TH=TH, method=method)
            L₁[n,1], L₂[n,1], SNR[n,1], PSNR[n,1] = getresults(X₀,Y₀,Y₁,TH,th)
            L₁[n,2], L₂[n,2], SNR[n,2], PSNR[n,2] = getresults(X₀,Y₀,Y₁,TH,mean(th))
            L₁[n,3], L₂[n,3], SNR[n,3], PSNR[n,3] = getresults(X₀,Y₀,Y₁,TH,median(th))
        end
    end

    return mean(L₁,dims=1), mean(L₂,dims=1), mean(SNR,dims=1), mean(PSNR,dims=1)
end

## Shifted experiment 
L1, L2, SNR, PSNR = L8testdenoise(testdata.blocks, 1.25, TH=SoftTH(), method=:jeff);

## Non-shifted experiment
L1, L2, SNR, PSNR = L8testdenoise_noshift(testdata.mishmash, 0.7, TH=SoftTH(), method=:jeff);

## Export Results
df = DataFrame(L1 = vec(L1), L2 = vec(L2), PSNR = vec(PSNR), SNR = vec(SNR))
CSV.write("../../Desktop/results.csv", df)
