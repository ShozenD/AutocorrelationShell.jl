module ACWT
export
    acfilter,
    leftchild,
    rightchild
using Wavelets, LinearAlgebra

"""Computes ac coefficients of given filter"""
function autocorr(f::OrthoFilter)
    H = WT.qmf(f)
    l = length(H)
    result = zeros(l - 1)
    for k in 1:(l - 1)
        for i in 1:(l - k)
            @inbounds result[k] += H[i] * H[i + k]
        end
        result[k] *= 2
    end
    return result
end

"""Computes low ac filter"""
function pfilter(f::OrthoFilter)
    a = autocorr(f)
    c1 = 1 / sqrt(2)
    c2 = c1 / 2
    b = c2 * a
    return vcat(b[end:-1:1], c1, b)
end

"""Computes the high AC filter"""
function qfilter(f::OrthoFilter)
    a = autocorr(f)
    c1 = 1 / sqrt(2)
    c2 = c1 / 2
    b = -c2 * a
    return vcat(b[end:-1:1], c1, b)
end

function makeqmfpair(f::OrthoFilter)
    pmfilter, qmfilter = pfilter(f), qfilter(f)
    return pmfilter, qmfilter
end

function makereverseqmfpair(f::OrthoFilter)
    pmf, qmf = makeqmfpair(f)
    return reverse(pmf), reverse(qmf)
end

"""For acwpt array methods"""
leftchild(i::Integer) = i<<1;
rightchild(i::Integer) = i<<1+1;

end # End module
