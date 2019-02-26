using DSP

function autocorr(H)
    """autocorr: Compute the autocorrelation of a filter with itself.

    Input:
        H: Filter
    Output:
        Autocorrelation coefficients of H
"""
    l = length(H)
    result = zeros(1, l - 1)
    for k in 1:(l - 1)
        for i in 1:(l - k)
            result[k] = 2 * H[i] * H[i + k]
        end
    end

    return result
end
