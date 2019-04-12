function iwt_ac(acwt)
    w = acwt[:, 1]
    n, m = size(acwt)
    for i = 2:m
        y = (y + acwt[:, i]) / sqrt(2)
    end
    return y'
end
