function cirperm(a)
    n = length(a)
    res = zeros(n, n)
    for i = 1:n
        res[i, :] = translate(a, 1 - i)
    end
    return res
end
