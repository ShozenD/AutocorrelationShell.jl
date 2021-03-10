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

"""
jbb_costs(X::AbstractArray{<:Number, 2}; et::CostFunction=NormCost())
Computes Norm Entropy/Loglp cost of each node in a binary tree in order to find the joint best basis
Input:
    - `X::AbstractArray`: array of size (n, L, N)
    - `et::CostFunction`: cost function for the selection of joint best basis, defaults to NormCost()
Output:
    - Vector of length 2^levels - 1, each entry corresponds to the cost in the node
"""
function jbb_costs(X::Array{Float64,3}; et::Wavelets.Entropy=NormEntropy())
    # For each level, compute mean, sum of squares, and variance
    (n, L, N) = size(X)
    EX = sum(X, dims = 3) ./ N              # calculate E(X)
    EX² = sum(X .^ 2, dims = 3) ./ N        # calculate E(X²)
    VarX = EX² - (EX) .^ 2                  # calculate Var(X)
    VarX = reshape(VarX, (n, L))
    σ = VarX .^ 0.5                         # calculate σ = √Var(X)
    @assert all(σ .>= 0)

    costs = Vector{Float64}(undef, 2^L - 1)

    i = 1   # iterates over the nodes for the costs variable
    for lvl in 1:L
        node_len = node_length(n, lvl-1)
        for node in 1:(2^(lvl-1))
            coef = σ[((node-1)*node_len+1):(node*node_len), lvl]
            costs[i] = node_cost(coef, et)
            i += 1
        end
    end
    return costs
end

function jbbfindcost(X::Array{Float64,3}, et::Wavelets.Entropy=NormEntropy())
    nx,ny,nz = size(X)
    costmat = Array{Float64,2}(undef,(ny,nz))
    for idx in 1:nz
        @inbounds costmat[:,idx] = findcost(X[:,:,idx],et)
    end
    return mean(costmat,dims=2) # Take the average of all vectors
end

function jbbtree(X::Array{Float64,3}, et::Wavelets.Entropy=NormEntropy())
    M = size(X,2)
    costvec = jbbfindcost(X,et)
    besttree = trues(M)

    for idx in (M>>1):-1:1
        if costvec[idx] <= (costvec[idx<<1] + costvec[(idx<<1)+1])/2
            besttree[subtree(idx,M)] .= false
        else
            costvec[idx] = (costvec[idx<<1] + costvec[(idx<<1)+1])/2
        end
    end
    return besttree
end

function JBBΣ(X::Array{Float64,2}, wt::OrthoFilter, n::Integer, noise::Float64;
              TH::Wavelets.Threshold.THType=HardTH(), L::Integer=maxtransformlevels(X[:,1]),
              method::Symbol=:jeff)
    # add noise
    X₁ = addnoise(X,n,noise)
    # acwpt decompositions
    Y₀ = decompositions(X₁,wt,L)
    # Find best basis trees for each signal
    y₁ = jbbtree(Y₀)
    # Find threshold
    if TH == SoftTH()
        th = bestthreshold(Y₀,y₁,TH=TH,method=method,elbows=3)
    else
        th = bestthreshold(Y₀,y₁,TH=TH,method=method,elbows=2)
    end

    return X₁, Y₀, y₁, th
end

## Load data
PATH = "/Users/shozendan/Documents/autocorrelation-shell/test/data/wavelet_test_256.csv";
testdata = CSV.read(PATH, DataFrame)

## Load Data
function jbbtestdenoise(x::Vector{Float64}, noise::Float64;
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
            X₁, Y₀, Y₁, th = JBBΣ(X₀, wt, n, noise, TH=TH, method=method)
            L₁[n,1], L₂[n,1], SNR[n,1], PSNR[n,1] = getresults(X₀,Y₀,Y₁,TH,th)
            L₁[n,2], L₂[n,2], SNR[n,2], PSNR[n,2] = getresults(X₀,Y₀,Y₁,TH,mean(th))
            L₁[n,3], L₂[n,3], SNR[n,3], PSNR[n,3] = getresults(X₀,Y₀,Y₁,TH,median(th))
        end
    end

    return mean(L₁,dims=1), mean(L₂,dims=1), mean(SNR,dims=1), mean(PSNR,dims=1)
end

function jbbtestdenoise_noshift(x::Vector{Float64}, noise::Float64;
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
            X₁, Y₀, Y₁, th = JBBΣ(X₀, wt, n, noise, TH=TH, method=method)
            L₁[n,1], L₂[n,1], SNR[n,1], PSNR[n,1] = getresults(X₀,Y₀,Y₁,TH,th)
            L₁[n,2], L₂[n,2], SNR[n,2], PSNR[n,2] = getresults(X₀,Y₀,Y₁,TH,mean(th))
            L₁[n,3], L₂[n,3], SNR[n,3], PSNR[n,3] = getresults(X₀,Y₀,Y₁,TH,median(th))
        end
    end

    return mean(L₁,dims=1), mean(L₂,dims=1), mean(SNR,dims=1), mean(PSNR,dims=1)
end

L1,L2,SNR,PSNR = jbbtestdenoise(testdata.blocks, 0.5, TH=SoftTH(), method=:jeff);

L1,L2,SNR,PSNR = jbbtestdenoise_noshift(testdata.blocks, 1.25, TH=HardTH(), method=:jeff);

L1

L2

PSNR

SNR
