function cdthresh(x, th)
    return diag(abs(x) >= th) * x
end
