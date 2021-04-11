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
            i = mod1(i + (L÷2+1) * 2^(j-1),N) # Need to shift by half the filter size because of periodicity assumption 
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
    Pmf, Qmf = ACWT.makeqmfpair(wt)
	wp = zeros(n,L+1)
	wp[:,1] = x

    for j in 1:L
        @inbounds wp[:,1], wp[:,L+2-j] = acwt_step(wp[:,1],j,Qmf,Pmf)
    end

	return wp
end

x = zeros(128)
x[64] = 1

acwt(x, wavelet(WT.haar)) |> wiggle

acwt(x, wavelet(WT.db3)) |> wiggle

acwt(x, wavelet(WT.db5)) |> wiggle

