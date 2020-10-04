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

"""
	dyadlength(x::Vector{<:Number})

Returns the dyadic length of a sequence `x`
"""
function dyadlength(x::Vector{<:Number})
    return trunc(Integer, log2(length(x)))
end

"""
	echant(n::Integer, d::Integer, b::Integer)

Sample the range `(b+1):n` in intervals of `2^d`
"""
function echant(n::Integer, d::Integer, b::Integer)
    return (b + 1):(2^d):n
end

"""
	autocorr(H::AbstractArray{<:Number})

Computes the autocorrelation coefficients of a given filter
"""
function autocorr(H::AbstractArray{<:Number})
    l = length(H)
    result = zeros(l - 1)
    @inbounds begin
        for k in 1:(l - 1)
            for i in 1:(l - k)
                result[k] += H[i] * H[i + k]
            end
            result[k] *= 2
        end
    end
    return result
end

"""
	autocorr(f::OrthoFilter)

Computes the autocorrelation coefficients of a given filter
"""
function autocorr(f::OrthoFilter)
    return autocorr(WT.qmf(f))
end

"""
    pfilter(filter::OrthoFilter)

Computes the low AC filter.

# Example
```julia
using Wavelets
H = wavelet(WT.db2)
P = pfilter(H)
```
"""
function pfilter(filter::OrthoFilter)
    a = autocorr(filter)
    c1 = 1 / sqrt(2)
    c2 = c1 / 2
    b = c2 * a
    return vcat(b[end:-1:1], c1, b)
end

"""
	qfilter(filter::OrthoFilter)

Computes the high AC filter.

# Example
```julia
using Wavelets
H = wavelet(WT.db2)
Q = qfilter(H)
```
"""
function qfilter(filter::OrthoFilter)
    a = autocorr(filter)
    c1 = 1 / sqrt(2)
    c2 = c1 / 2
    b = -c2 * a
    return vcat(b[end:-1:1], c1, b)
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

The inverse function of `fwt_ac`. Reconstructs the original signal given an array of autocorrelation wavelet coefficients.
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

    # Sanity Check
    @assert L >= 0
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

"""
    acwt(x; L, P, Q)

Computes the forward autocorrelation wavelet transform. Wrapper for the fwt_ac function.

# Arguments
- `x::Vector{<:Real}`: array to transform.
- `L::Integer`: degree of coarsest scale. *default* = 1.
- `P::Vector{<:Real}`: Low AC shell filter.
- `Q::Vector{<:Real}`: High AC shell filter.

# Example
```julia
using Wavelets

H = wavelet(WT.db2)
Q = qfilter(H)
P = pfilter(H)

x = zeros(256)
x[128] = 1

decomp = acwt(x; P=P, Q=Q)
```
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
