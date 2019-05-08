module AutocorrelationShell

using DSP
using StatsBase
using Wavelets

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
   ypadded = filter(f, xpadded)
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

function dyadlength(x)
"""

"""
    return log2(length(x))
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
    result = zeros(1, l - 1)
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

function Pfilter(filter::OrthoFilter)
"""
	Pfilter(filter::OrthoFilter)

	Computes the autocorrelation low filter
"""
    a = autocorr(filter)
    c1 = 1 / (2 * sqrt(2))
    c2 = c1 / 2
    b = c2 * a
    return vcat(b[end:-1:1], c1, b)
end

function Qfilter(filter::OrthoFilter)
"""
	Qfilter(filter::OrthoFilter)

	Computes the autocorrelation high filter
"""
    a = autocorr(filter)
    c1 = 1 / (2 * sqrt(2))
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
    p = length(Q)
    tran2 = p - 1
    tran1 = tran2 ÷ 2

    d = translate(iconv(Q, translate(x, -tran1)), tran2)
    return d[1:n]
end

function iwt_ac(acwt)
"""
	iwt_ac(acwt)

	Performs the 1D autocorrelation decomposition (inverse)
"""
    w = acwt[:, 1]
    n, m = size(acwt)
    for i = 2:m
        y = (y + acwt[:, i]) / sqrt(2)
    end
    return y'
end

function fwt_ac(x,L,P,Q)
"""
	fwt_ac(x,L,P,Q)

	Computes the forward autocorrelation wavelet transform

	# Arguments
	- `x`: array to transform
	- `L`: degree of coarsest scale
	- `P`: Low AC shell filter
	- `Q`: High AC shell filter
"""
	n = length(x)
	J = dyadlength(n)
	D = J-L
	wp = zeros(n,D+1)
	x = reshape(x,(1,length(x)))

	wp[:,1] = x'
	for d=0:(D-1)
		for b=0:(2^d-1)
		   s = wp[echant(n,d,b),1]'
		   h = ac_filter(s,Q)
		   l = ac_filter(s,P)
		   wp[echant(n,d,b),D+1-d] = h'
		   wp[echant(n,d,b),1] = l'
		 end
	end
	return wp
end

function thresh0(acwt, th, hard=false)
"""
	thresh0(acwt, th, hard=false)

	Thresholds a computed signal decomposition (`acwt`) to threshold `th`. Either hard or soft (default) thresholding can be used
"""
	n, m = size(acwt)
	res=zeros(n,m)
	res[:,1]=acwt[:,1]
	if hard
		for i=2:m
			res[:,i] = threshold(acwt[:,i],HardTH(),th)
		end
	else
		for i=2:m
			res[:,i] = threshold(acwt[:,i],SoftTH(),th)
		end
	end
	return res
end

function ACCorrCalc(R, w::OrthoFilter, L)
"""
	ACCorrCalc(R, w::OrthoFilter, L)

	Calculates all induced autocovariance functions at all levels

	# Arguments
	- `R`: autocovariance of the original signal
	- `w`: wavelet to be used
	- `L`: decomposition level
"""
  P = Pfilter(w)
  Q = Qfilter(w)

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

    return tab2[:, 1]'
end


export autocorr

end # module
