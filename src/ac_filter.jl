function ac_filter(x, filter)
    n = length(x)
    p = length(Q)
    tran2 = p - 1
    tran1 = tran2 ÷ 2

    d = translate(iconv(Q, translate(x, -tran1)), tran2)
    return d[1:n]
end
