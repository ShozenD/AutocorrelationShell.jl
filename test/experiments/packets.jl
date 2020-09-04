include("../src/AutocorrelationShell.jl")
using Wavelets, Main.AutocorrelationShell

function dyadlength(x)
"""
    dyadlength(x)

    Returns the dyadic length of a sequence `x`
"""
    return trunc(Interger, log2(length(x)))
end

function binary_ac1d(x, L, P, Q, wp_array, i, depth, max_depth)
    if depth <= max_depth
        wp_array[:, i] = x
        decomp = fwt_ac(x, L, P, Q)
        
        # left node
        left = decomp[:,1]
        binary_ac1d(left, L, P, Q, wp_array, 2*i, depth+1, max_depth)
        
        # right node
        right = decomp[:,2]
        binary_ac1d(right, L, P, Q, wp_array, 2*i+1, depth+1, max_depth)
    end
end

function ac1d_wpt(x, P, Q)
    max_depth = dyadlength(x)
    max_node_index = 2^max_depth - 1
    L = max_depth - 1
    wpt_decomp = Array{Float64, 2}(undef, length(x), max_node_index)
    binary_ac1d(x, L, P, Q, wpt_decomp, 1, 1, max_depth)
    
    return wpt_decomp
end

x = rand(16)
H = wavelet(WT.db2)
Q = qfilter(H)
P = pfilter(H)
wpt_decomp = ac1d_wpt(x, P, Q)