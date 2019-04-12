using Wavelets

function Pfilter(filter::OrthoFilter)
    a = autocorr(filter)
    c1 = 1 / (2 * sqrt(2))
    c2 = c1 / 2
    b = c2 * a
    return vcat(b[end:-1:1], c1, b)
end
