function inv_ac_table(table, basis)
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
