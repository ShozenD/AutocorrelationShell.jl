using DSP

function autocorr(H::AbstractArray)
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
            result[k] += H[i] * H[i + k]
        end
        result[k] *= 2
    end

    return result
end

function autocorr(f::OrthoFilter)
    return autocorr(WT.qmf(f))
end
