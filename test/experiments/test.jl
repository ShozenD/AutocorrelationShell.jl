using AutocorrelationShell, LinearAlgebra, Wavelets, Plots, BenchmarkTools

function acwt_step(v::AbstractVector{T}, j::Integer, h::Array{S,1},
        g::Array{S,1}) where {T <: Number, S <: Number}
    N = length(v)
    L = length(h)
    v1 = zeros(T, N)
    w1 = zeros(T, N)
    @inbounds begin
        for i in 1:N
            t = i
            if N ÷ 2^(j-1) > L # adjusting index to deal with shift
                i = mod1(i+2^(j+2),N)
            end
            for n in 1:L
                t += 2^(j-1)
                t = mod1(t, N)
                w1[i] += h[n] * v[t]
                v1[i] += g[n] * v[t]
            end
        end
    end
    return v1, w1
end

function acwt_test(x::AbstractVector{<:Number}, wt::OrthoFilter, L::Integer=maxtransformlevels(x))

    @assert L <= maxtransformlevels(x) || throw(ArgumentError("Too many transform levels (length(x) < 2^L"))
    @assert L >= 1 || throw(ArgumentError("L must be >= 1"))

    # Setup
	n = length(x)
    Pmf, Qmf = ACWT.makereverseqmfpair(wt)
	wp = zeros(n,L+1)
	wp[:,1] = x

    for j in 1:L
        @inbounds wp[:,1], wp[:,L+2-j] = acwt_step(wp[:,1],j,Qmf,Pmf)
    end

	return wp
end

# Array method
function acwpt_test(W::Array{<:Number,2}, i::Integer, d::Integer,
               Pmf::Vector{<:Number}, Qmf::Vector{<:Number})
    n,m = size(W)
    if i<<1+1 <= m
        W[:,i<<1], W[:,i<<1+1] = acwt(W[:,i],d,Qmf,Pmf)
        acwpt_test(W,i<<1,d+1,Pmf,Qmf) # left
        acwpt_test(W,i<<1+1,d+1,Pmf,Qmf) # right
    end
end

# Combined
function acwpt_test(x::Vector{T}, wt::OrthoFilter,
               L::Integer=maxtransformlevels(x), method::Symbol=:array) where T<:Number
    if method == :tree
        root = AcwptNode(x)
        acwpt(x,root,wt,L)
        return root
    elseif method == :array
        W = Array{Float64,2}(undef,length(x),2^(L+1)-1)
        W[:,1] = x
        Pmf, Qmf = ACWT.makereverseqmfpair(wt)
        acwpt_test(W,1,1,Pmf,Qmf)
        return W
    else
        throw(ArgumentError("unkown method"))
    end
end

x = randn(128)
x[64] = 20

y1 = acwpt_test(x, wavelet(WT.db4))

y2 = acwpt(x, wavelet(WT.db4), maxtransformlevels(x), :array)

norm(y1[:,128]-y2[:,128])
