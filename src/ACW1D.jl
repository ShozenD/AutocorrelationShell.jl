function conv_filt(a,b,sig)
    rst = zeros(size(sig))
    for i in eachindex(sig)
        tmp =  sum(b[j] * sig[i-j+1] for j in 1:min(i, length(b)))
        tmp -= sum(a[j] * rst[i-j+1] for j in 1:min(i, length(a)))
        rst[i] = tmp / a[1]
    end
    return rst
end

function iconv(f,x)
"""

"""
   n = length(x)
   p = length(f)
   if p <= n
      xpadded = vcat(x[(n+1-p):n], x)
   else
      z = zeros(p)
      for i=1:p
          imod = 1 + rem(p*n -p + i-1,n)
          z[i] = x[imod]
      end
      xpadded = vcat(z,x)
   end
   ypadded = conv_filt(1, f, xpadded)
   y = ypadded[(p+1):(n+p)]
   return y
end

function translate(x, to)
"""
	translate(x, to)

	Circular translates the signal `x` by the desired places `to` to the left
"""
    return circshift(x, -to)
end

function dyadlength(x::Vector{<:Number})
"""
	dyadlength(x)

	Returns the dyadic length of a sequence `x`
"""
    return trunc(Integer, log2(length(x)))
end

function subsample(x)
"""
	subsample(x)

	Samples every other element of input signal `x`
"""
    return x[1:2:end-1]
end

function cdthresh(x, th)
"""
	cdthresh(x, th)

	Thresholds the elements of `x` by threshold `th`
"""
    return diag(abs(x) >= th) * x
end

function cirperm(a)
"""
	cirperm(a)

	Given a sequence `a` for length `n`, return an `n`x`n` matrix containing all circular permutations of `a`
"""
    n = length(a)
    res = zeros(n, n)
    for i = 1:n
        res[i, :] = translate(a, 1 - i)
    end
    return res
end

function node(d, b)
"""

"""
  return 2^d + b
end

function echant(n, d, b)
"""
	echant(n, d, b)

	Sample the range `(b+1):n` in intervals of `2^d`
"""
    return (b + 1):(2^d):n
end

function autocorr(H::AbstractArray)
"""
	autocorr(H)

	Computes the autocorrelation coefficients of a given filter
"""
    l = length(H)
    result = zeros(l - 1)
    for k in 1:(l - 1)
        for i in 1:(l - k)
            result[k] += H[i] * H[i + k]
        end
        result[k] *= 2
    end

    return result
end

function autocorr(f::OrthoFilter)
"""
	autocorr(f)

	Computes the autocorrelation coefficients of a given filter
"""
    return autocorr(WT.qmf(f))
end

function pfilter(filter::OrthoFilter)
"""
	Pfilter(filter::OrthoFilter)

	Computes the autocorrelation low filter
"""
    a = autocorr(filter)
    c1 = 1 / sqrt(2)
    c2 = c1 / 2
    b = c2 * a
    return vcat(b[end:-1:1], c1, b)
end

function qfilter(filter::OrthoFilter)
"""
	Qfilter(filter::OrthoFilter)

	Computes the autocorrelation high filter
"""
    a = autocorr(filter)
    c1 = 1 / sqrt(2)
    c2 = c1 / 2
    b = -c2 * a
    return vcat(b[end:-1:1], c1, b)
end

function acnyquist(s)
"""
	acnyquist(s)

	Computes the Nyquist frequency of a given signal `s`
"""
	n = length(s)
	sub = s[collect(2:2:n)]
	r = 0
	n = length(sub)
	for k = 1:n
		r += sub[k]*((-1)^(k-1))
	end
	return r
end

function ac_filter(x, filter)
"""
	ac_filter(x, filter)

	Computes the response of signal `x` to the autocorrelation filter `filter`
"""
    n = length(x)
    p = length(filter)
    tran2 = p - 1
    tran1 = tran2 ÷ 2

    d = circshift(iconv(filter, circshift(x, tran1)), -tran2)
    return d[1:n]
end

function iwt_ac(acwt::AbstractArray{<:Number})
"""
	iwt_ac(acwt)

Performs the 1D autocorrelation decomposition (inverse)
"""
    y = deepcopy(acwt[:, 1])
    n, m = size(acwt)
    @inbounds begin
        for i = 2:m
            y = (y + acwt[:, i]) / sqrt(2)
        end
    end
    return y
end

"""
	fwt_ac(x, L, P, Q)

Computes the forward autocorrelation wavelet transform for a given signal.

# Arguments
- `x::Vector{<:Number}`: Signal.
- `L::Integer`: Degree of coarsest scale.
- `P::Vector{<:Number}`: Low AC shell filter.
- `Q::Vector{<:Number}`: High AC shell filter.
"""
function fwt_ac(x::Vector{T}, L::Integer, P::Vector{T}, Q::Vector{T}) where T<:Number

	n = length(x)
	J = dyadlength(x)

    @assert L > 0
    @assert L <= J

	D = J-L
	wp = zeros(n,D+1)

	wp[:,1] = x
    @inbounds begin
    	for d=0:(D-1)
    		for b=0:(2^d-1)
    		   s = wp[echant(n,d,b),1]
    		   h = ac_filter(s,Q)
    		   l = ac_filter(s,P)
    		   wp[echant(n,d,b),D+1-d] = h
    		   wp[echant(n,d,b),1] = l
    		 end
    	end
    end
	return wp
end

function autocorr_calc(R, w::OrthoFilter, L)
"""
	autocorr_calc(R, w::OrthoFilter, L)

	Calculates all induced autocovariance functions at all levels

	# Arguments
	- `R`: autocovariance of the original signal
	- `w`: wavelet to be used
	- `L`: decomposition level
"""
	P = pfilter(w)
	Q = qfilter(w)

	b = autocorr(P) / 2
	c = autocorr(Q) / 2

	b = vcat(b[end:-1:1], norm(P)^2, b)
	c = vcat(c[end:-1:1], norm(Q)^2, c)

	n = length(R)
	J = log2(n)

	RS = zeros(n, L+1)
	RD = zeros(n, L)
	RS[:,1] = R

	for j = 1:L
	RS[1:2^(J - j), j + 1] = subsample(ac_filter(RS[1:2^(J - j + 1), j]', b))';
	RD[1:2^(J - j), j] = subsample(ac_filter(RS[1:2^(J - j + 1), j]',c))';
	end

	return RS, RD
end

function inv_ac_table(table, basis)
"""
	inv_ac_table(table, basis)

	Computes the inverse AC decomposition using the chosen basis
"""
    n, D = size(table)
    L = floor(log2(D))
    tab2 = deepcopy(table)

    for d = (L-1):-1:0
        for b = 0:(2^d - 1)
            if basis[node(d, b)] == 1
                tab2[:, node(d, b)] = iwt_ac(tab2[:, node(d + 1, 2 * b):node(d + 1, 2 * b + 1)])'
			end
        end
    end

    return tab2[:, 1]'
end

"""
    acwt(x, L, P, Q)

Computes the forward autocorrelation wavelet transform. Wrapper for the fwt_ac function.

# Arguments
- `x::Vector{<:Real}`: array to transform.
- `L::Integer`: degree of coarsest scale. *default*: 1.
- `P::Vector{<:Real}`: Low AC shell filter.
- `Q::Vector{<:Real}`: High AC shell filter.
"""
function acwt(x::Vector{T}; L::Integer=1, P::Vector{T}, Q::Vector{T}) where T<:Real
    return fwt_ac(x, L, P, Q)
end

"""
    iacwt(acwt)

Inverse autocorrelation wavelet transform(signal reconstruction). Wrapper for the iwt_ac function.

# Arguments
- `acwt::AbstractArray{<:Number}`: Array of wavelet coefficients.
"""
function iacwt(acwt::AbstractArray{<:Number})
    return iwt_ac(acwt)
end
