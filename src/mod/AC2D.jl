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
        accoef_matrix_3d[:,:,i] = fwt_ac(img[:,i],L,P,Q)
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
    accoef_matrix_3d = ac1d_col(img, L, P, Q)
    accoef_matrix_4d = Array{Float64, 4}(undef, D, num_row, D, num_col)
    for i in 1:D
        accoef_matrix_4d[i,:,:,:] = ac1d_col(accoef_matrix_3d[:,i,:],L,P,Q)
    end
    accoef_matrix_4d = permutedims(accoef_matrix_4d, [2,4,1,3])
    return accoef_matrix_4d
end

# ------------ Inverse functions -------------
function rowinv(x)
"""
    rowinv(x)
Performs iwt_ac on rows of each element of the input
vector.

# Arguments
- `x`: one element of the ac2d function output vector
"""
    n_scale=length(x)
    num_row,num_col=size(x[1])
    a = [Array{Float64}(undef,num_row,n_scale) for i=1:num_col]
    for i=1:num_row
        for j=1:n_scale
            a[i][:,j]=x[j][i,:]
        end
    end
    b=[iwt_ac(a[i]) for i=1:length(a)]
    return transpose(reduce(hcat,b))
end



function iac2d(ac)
"""
    iac2d(ac)

Performs the inverse 2D autocorrelation wavelet transform.

# Arguments
- `ac`: ac2d function output
"""
    tmp = [rowinv(ac[i]) for i=1:length(ac)]
    n_scale=length(tmp)
    num_row,num_col=size(tmp[1])
    a = [Array{Float64}(undef,num_row,n_scale) for i=1:num_col]
    for i=1:num_row
        for j=1:n_scale
            a[i][:,j]=tmp[j][:,i]
        end
    end
    b=[iwt_ac(a[i]) for i=1:length(a)]
    return reduce(hcat,b)
end

end # module
