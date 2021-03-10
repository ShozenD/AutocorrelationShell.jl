function recurve(coef::AbstractArray{Float64,2}, bt::BitArray{1}, TH::Wavelets.Threshold.THType)
    f = iacwpt(coef,bt) # Original signal
    leaves = collectleaves(bt,true)
    leafnum = length(leaves) # number of leaf nodes
    if leafnum == 0 # No interesting decomposition found
        return 0, 0 # Relative error is always zero
    else
        α = coef[:,leaves]
        β = coef[:,leaves]
        maxcoef = maximum(abs.(α))
        τ = [0:0.01:round(maxcoef,digits=2);]
        ϵ = Vector{Float64}(undef,length(τ))
        for (idx,th) in enumerate(τ)
            coef[:,leaves] = threshold!(β,TH,th)
            g = iacwpt(coef,bt)
            ϵ[idx] = relerror(f,g)
            coef[:,leaves] = α
        end
        return τ, ϵ
    end
end

relerror(f::Vector{Float64},g::Vector{Float64},p::Integer=2) = norm(f-g,p)/norm(f,p)

## THRESHOLD SELECTION
# Determine best threshold for one signal
function bestthreshold(coef::AbstractArray{T}, bt::BitArray{1};
                       TH::Wavelets.Threshold.THType=HardTH(), elbows::Integer=2,
                       method::Symbol=:shozen, plot::Bool=false) where T<:Float64
    @assert elbows >= 1
    @assert method ∈ [:shozen, :jeff]

    if method === :shozen
        x,r = recurve(coef, bt, TH)
        if x == 0 return 0 end # No interesting decomposition found
        x = sort(abs.(x), rev=true)
        r = sort(r, rev=true)
    else
        leaves = collectleaves(bt,true)
        if length(leaves) == 0 # No interesting decomposition found
            α = coef[:,1] # use the input signal
        else
            α = collect(Iterators.flatten(coef[:,leaves]))
        end
        x = sort(abs.(α), rev = true)
        r = orth2relerror(α)
    end

    # shift the data points
    push!(x, 0)
    pushfirst!(r, r[1])

    # reorder the data points
    xmax = maximum(x)
    ymax = maximum(r)
    x = x[end:-1:1]/xmax
    y = r[end:-1:1]/ymax

    # compute elbow method
    ix = Vector{Int64}(undef, elbows)
    A = Vector{Vector{Float64}}(undef, elbows)
    v = Vector{Vector{Float64}}(undef, elbows)
    for i in 1:elbows
        (ix[i], A[i], v[i]) = i > 1 ? findelbow(x[1:ix[i-1]], y[1:ix[i-1]]) : findelbow(x, y)
    end

    if plot == true
        p = threshold_vs_relativeerror_plot(x*xmax, y*ymax, ix, A, v)
        display(p)
    end

    return x[ix[end]]*xmax
end

"""
orth2relerror(orth::Vector{T}) where T<:Number
Given a vector 'orth' of orthonormal expansion coefficients, return a
vector of relative approximation errors when retaining the 1,2,...,N
largest coefficients in magnitude.
Input
   orth        a vector of orthonormal expansion coefficients
Output
   relerror    a vector of relative approximation errors
"""
function orth2relerror(orth::AbstractVector{T}) where T<:Number
    # sort the coefficients
    orth = sort(orth.^2, rev = true)

    # compute the relative errors
    return ((abs.(sum(orth) .- cumsum(orth))).^0.5) / sum(orth).^0.5
end

# Determine best threshold for each signal
function bestthreshold(coef::Array{Float64,3}, BT::BitArray{2};
                              TH::Wavelets.Threshold.THType=HardTH(),
                              elbows=2, method::Symbol=:shozen,
                              plot::Bool=false)
    thres = Vector{Float64}(undef,size(coef,3))
    for idx in axes(coef,3)
        @inbounds thres[idx] = bestthreshold(coef[:,:,idx],BT[:,idx],TH=TH,elbows=elbows, method=method)
    end

    return thres
end

function bestthreshold(coef::Array{Float64,3}, BT::BitArray{1};
                              TH::Wavelets.Threshold.THType=HardTH(),
                              elbows=2, method::Symbol=:shozen,
                              plot::Bool=false)
    thres = Vector{Float64}(undef,size(coef,3))
    for idx in axes(coef,3)
        @inbounds thres[idx] = bestthreshold(coef[:,:,idx],BT,TH=TH,elbows=elbows,method=method)
    end

    return thres
end

"""
findelbow(x::Vector{T}, y::Vector{T}, index::Integer) where T<:Number
Given the x and y coordinates of a curve, find the elbow.
Input:
   - x       the x coordinates of the points
   - y       the y coordinates of the points
   - index   the index of the endpoint
Output:
   - IX      the index of the elbow
"""
function findelbow(x::AbstractVector{T}, y::AbstractVector{T}) where T<:Number
    # a unit vector pointing from (x1,y1) to (xN,yN)
    v = [x[end] - x[1], y[end] - y[1]]
    v = v/norm(v,2)

    # subtract (x1,y1) from the coordinates
    xy = [x.-x[1] y.-y[1]]

    # the hypothenuse
    H = reshape((sum(xy.^2, dims = 2)).^0.5, :)

    # the adjacent side
    A = xy * v

    # the opposite side
    O = abs.(H.^2 - A.^2).^0.5       # TODO: fix case where value approaches machine accuracy

    # return the largest distance
    return (findmax(O)[2], A, v)
end

# threshold vs relative error curve
function threshold_vs_relativeerror_plot(x::Vector{<:Number}, y::Vector{<:Number}, ix::Vector{<:Integer}, A::Vector{<:Vector{<:Number}}, v::Vector{<:Vector{<:Number}})
    # rescale x and y values
    xmax = maximum(x)
    ymax = maximum(y)

    # relative error line
    p = plot(x, y, lw = 2, color = :blue, legend = false)
    plot!(p, xlims = (0, 1.004*xmax), ylims = (0, 1.004*ymax))

    @assert length(ix) == length(A) == length(v)
    elbows = length(ix)
    for i in 1:elbows
        color = i + 1
        # diagonal line
        endpoint = i > 1 ? ix[i-1] : length(x)
        plot!([x[1], 1.004*x[endpoint]], [y[1], 1.004*y[endpoint]], lw = 2, color = color)
        # perpendicular line
        dropto = [x[1], y[1]] + A[i][ix[i]]*(v[i].*[xmax, ymax])
        plot!(p, [x[ix[i]], dropto[1]], [y[ix[i]], dropto[2]], lw = 2, color = color)
        # highlight point
        scatter!(p, [x[ix[i]]], [y[ix[i]]], color = color)
    end

    # add plot labels
    plot!(p, xlabel = "Threshold", ylabel = "Relative Error")
    return p
end


## Helper functions
function psnr(x::AbstractVector{T},x₀::AbstractVector{T}) where T<:Number
    sse = zero(T)
    for i in eachindex(x)
        @inbounds sse += (x[i] - x₀[i])^2
    end
    mse = sse/length(x)
    return 20 * log(10, maximum(x₀)) - 10 * log(10, mse)
end

function noshiftsignal(x::Vector{Float64}, M::Integer)
    y = Array{Float64,2}(undef, length(x), M)
    for idx in 1:100
        @inbounds y[:,idx] = x
    end
    return y
end

function shiftsignal(x::Vector{Float64}, M::Integer)
    y = Array{Float64,2}(undef, (length(x), M)) # Circshifted signals
    for idx in 0:(M-1)
        @inbounds y[:,idx+1] = circshift(x,2*idx)
    end
    return y
end

function addnoise(X::Array{Float64,2}, n::Integer, a::Float64)
    rng = MersenneTwister(n^2)
    return makenoisy(X, rng, a)
end

function decompositions(X::Array{Float64,2},wt::OrthoFilter,L::Integer)
    nx, nz = size(X)
    ny = nx<<1-1
    y = Array{Float64,3}(undef,(nx,ny,nz))
    for idx in 1:nz
        @inbounds y[:,:,idx] = acwpt(X[:,idx],wt,L,:array)
    end
    return y
end

function besttrees(X::Array{Float64,3}, et::Wavelets.Entropy=NormEntropy())
    nx,ny,nz = size(X)
    y = BitArray{2}(undef,(ny,nz))
    for idx in 1:nz
        @inbounds y[:,idx] = bestbasistree(X[:,:,idx],et)
    end
    return y
end

function Σ(X::Array{Float64,2}, wt::OrthoFilter, n::Integer, noise::Float64;
           TH::Wavelets.Threshold.THType=HardTH(), L::Integer=maxtransformlevels(X[:,1]),
           method::Symbol=:individual)
    # add noise
    X₁ = addnoise(X,n,noise)
    # acwpt decompositions
    Y₀ = decompositions(X₁,wt,L)
    # Find best basis trees for each signal
    Y₁ = besttrees(Y₀)
    # Find threshold
    if TH == SoftTH()
        th = bestthreshold(Y₀,Y₁,TH=TH,method=method,elbows=3)
    else
        th = bestthreshold(Y₀,Y₁,TH=TH,method=method,elbows=2)
    end

    return X₁, Y₀, Y₁, th
end

function getresults(X::Array{T,2}, Y::Array{T,3}, BT::BitArray{2},
                    TH::Wavelets.Threshold.THType, th::Array{T,1}) where T<:Float64
    M = size(Y,3)
    L₁ = Vector{Float64}(undef,M)
    L₂ = Vector{Float64}(undef,M)
    SNR = Vector{Float64}(undef,M)
    PSNR = Vector{Float64}(undef,M)

    for idx in 1:M
        leaves = collectleaves(BT[:,idx],true)
        coef = Y[:,:,idx]
        α = coef[:,leaves]
        coef[:,leaves] = threshold!(α,TH,th[idx])
        g = iacwpt(coef,BT[:,idx])
        L₁[idx] = relerror(X[:,idx],g,1)
        L₂[idx] = relerror(X[:,idx],g,2)
        SNR[idx] = snr(X[:,idx],g)
        PSNR[idx] = psnr(X[:,idx],g)
    end

    return mean(L₁), mean(L₂), mean(SNR), mean(PSNR)
end

function getresults(X::Array{T,2}, Y::Array{T,3}, BT::BitArray{2},
                    TH::Wavelets.Threshold.THType, th::Float64) where T<:Float64
    M = size(Y,3)
    L₁ = Vector{Float64}(undef,M)
    L₂ = Vector{Float64}(undef,M)
    SNR = Vector{Float64}(undef,M)
    PSNR = Vector{Float64}(undef,M)

    for idx in 1:M
        leaves = collectleaves(BT[:,idx],true)
        coef = Y[:,:,idx]
        α = coef[:,leaves]
        coef[:,leaves] = threshold!(α,TH,th)
        g = iacwpt(coef,BT[:,idx])
        L₁[idx] = relerror(X[:,idx],g,1)
        L₂[idx] = relerror(X[:,idx],g,2)
        SNR[idx] = snr(X[:,idx],g)
        PSNR[idx] = psnr(X[:,idx],g)
    end

    return mean(L₁), mean(L₂), mean(SNR), mean(PSNR)
end

function getresults(X::Array{T,2}, Y::Array{T,3}, BT::BitArray{1},
                    TH::Wavelets.Threshold.THType, th::Array{T,1}) where T<:Float64
    M = size(Y,3)
    L₁ = Vector{Float64}(undef,M)
    L₂ = Vector{Float64}(undef,M)
    SNR = Vector{Float64}(undef,M)
    PSNR = Vector{Float64}(undef,M)

    for idx in 1:M
        leaves = collectleaves(BT,true)
        coef = Y[:,:,idx]
        α = coef[:,leaves]
        coef[:,leaves] = threshold!(α,TH,th[idx])
        g = iacwpt(coef,BT)
        L₁[idx] = relerror(X[:,idx],g,1)
        L₂[idx] = relerror(X[:,idx],g,2)
        SNR[idx] = snr(X[:,idx],g)
        PSNR[idx] = psnr(X[:,idx],g)
    end

    return mean(L₁), mean(L₂), mean(SNR), mean(PSNR)
end

function getresults(X::Array{T,2}, Y::Array{T,3}, BT::BitArray{1},
                    TH::Wavelets.Threshold.THType, th::Float64) where T<:Float64
    M = size(Y,3)
    L₁ = Vector{Float64}(undef,M)
    L₂ = Vector{Float64}(undef,M)
    SNR = Vector{Float64}(undef,M)
    PSNR = Vector{Float64}(undef,M)

    for idx in 1:M
        leaves = collectleaves(BT,true)
        coef = Y[:,:,idx]
        α = coef[:,leaves]
        coef[:,leaves] = threshold!(α,TH,th)
        g = iacwpt(coef,BT)
        L₁[idx] = relerror(X[:,idx],g,1)
        L₂[idx] = relerror(X[:,idx],g,2)
        SNR[idx] = snr(X[:,idx],g)
        PSNR[idx] = psnr(X[:,idx],g)
    end

    return mean(L₁), mean(L₂), mean(SNR), mean(PSNR)
end
