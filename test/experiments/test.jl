using AutocorrelationShell, LinearAlgebra, Wavelets, Plots

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

x = zeros(128)
x[64] = 1

acwt(x, wavelet(WT.db4)) |> wiggle

acwt_test(x, wavelet(WT.db4)) |> wiggle

acwt_test(x, wavelet(WT.haar)) |> wiggle


ACWT.makeqmfpair(wavelet(WT.db4))

Pmf_haar, Qmf_haar = ACWT.makeqmfpair(wavelet(WT.haar))

Pmf_db2, Qmf_db2 = ACWT.makeqmfpair(wavelet(WT.db2))

L, H = acwt_step(x, 1, Qmf_haar, Pmf_haar)

LL, LH = acwt_step(L, 2, Qmf_haar, Pmf_haar)

LLL, LLH = acwt_step(L, 3, Qmf_haar, Pmf_haar)


plot(LL)

argmax(L) # under shifts by 2

argmax(LL) # under shifts by 6 (relative undershift is 4?)

argmax(LLL) # u.s. by 10 (relative undershift is 6)
