module AC2D
export
    # 2D forward autocorrelation wavelet transform
    ac2d,
    # inverse 2D autocorrelation wavelet transform
    iac2d
using ..AC1D

"""
  ac1d_col(x, L, P, Q)

Computes the column wise 1D fwt_ac coeficients for 2D signals.

### Arguments
- `x`: 2D signal (n ,m) matrix
- `L`: Decomposition level
- `P`: Low AC shell filter
- `Q`: High AC shell filter

Returns a multidimensional array of (num_col, levels_of_decomposition)
"""
function ac1d_col(x, L, P, Q)
    num_row, num_col = size(x)
    J = trunc(Integer, log2(num_col))
    D = J - L + 1
    accoef_matrix_3d = Array{Float64, 3}(undef, num_col, D, num_col)
    for i in 1:num_col
        accoef_matrix_3d[:,:,i] = fwt_ac(x[:,i],L,P,Q)
    end
    return accoef_matrix_3d
end

# ------- Two dimensional functions --------
"""
  ac2d(x, L, P, Q)

Computes autocorrelation wavelet coeficients for 2D signals.

### Arguments
- `x`: 2D signals (n ,m) matrix
- `L`: Decomposition level
- `P`: Low AC shell filter
- `Q`: High AC shell filter

Returns a the multidimensional matrix of
(num_rows, num_cols, levels_of_decomp, levels_of_decomp) that stores
the coefficients of the decomposed image.
"""
function ac2d(x, L, P, Q)
    num_row, num_col = size(x)
    J = trunc(Integer, log2(num_col))
    D = J - L + 1
    accoef_matrix_3d = ac1d_col(x, L, P, Q)
    accoef_matrix_4d = Array{Float64, 4}(undef, D, num_row, D, num_col)
    for i in 1:D
        accoef_matrix_4d[i,:,:,:] = ac1d_col(accoef_matrix_3d[:,i,:],L,P,Q)
    end
    accoef_matrix_4d = permutedims(accoef_matrix_4d, [4,2,1,3])
    return accoef_matrix_4d
end

# ------------ Inverse functions -------------
"""
    iac2d(x)

Performs the inverse 2D autocorrelation wavelet transform.

# Arguments
- `x`: ac2d function output
"""
function iac2d(x)
    num_row, num_col, D_row, D_col = size(x)
    accoef_matrix_4d = permutedims(x, [3,2,4,1])
    accoef_matrix_3d = Array{Float64, 3}(undef, num_col, D_col, num_row)
    for i in 1:D_row
        for j in 1:num_row
            accoef_matrix_3d[:,i,j] = iwt_ac(accoef_matrix_4d[i,:,:,j])
        end
    end
    reconst = Array{Float64, 2}(undef, 512, 512)
    for i in 1:num_col
        reconst[:,i] = iwt_ac(accoef_matrix_3d[:,:,i])
    end
    return reconst
end

end # module
