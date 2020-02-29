module AC2D
export
    # 2D forward autocorrelation wavelet transform
    ac2d
using ..AC1D

"""
  ac1d_col(x, L, P, Q)

Computes the column wise 1D fwt_ac coeficients for 2D signals.

# Arguments
- `x`: 2D signal (n ,m) matrix
- `L`: Decomposition level
- `P`: Low AC shell filter
- `Q`: High AC shell filter
"""
function ac1d_col(x,L,P,Q)
  num_row, num_col = size(x)
  accoef_arr = [fwt_ac(x[:,i],L,P,Q) for i=1:num_col]
  n_scale = size(accoef_arr[1])[2]
  a = [Array{Float64}(undef,num_row,num_col) for i=1:n_scale]
  for i = 1:n_scale
      for j = 1:length(accoef_arr)
          a[i][:,j] = accoef_arr[j][:,i]
      end
  end
  return a
end

"""
  ac1d_row(x, L, P, Q)

Computes the row wise 1D fwt_ac coeficients for 2D signals.

# Arguments
- `x`: 2D signal (n ,m) matrix
- `L`: Decomposition level
- `P`: Low AC shell filter
- `Q`: High AC shell filter
"""
function ac1d_row(x,L,P,Q)
   xt = transpose(x)
   num_row,num_col=size(xt)
   accoef_arr = [fwt_ac(xt[:,i],L,P,Q) for i=1:num_col]
   n_scale = size(accoef_arr[1])[2]
   a = [Array{Float64}(undef,num_row,num_col) for i=1:n_scale]
   for i in 1:n_scale
       for j in 1:length(accoef_arr)
           a[i][j,:] = accoef_arr[j][:,i]
       end
   end
   return a
end

# ------- Two dimensional functions --------
"""
  ac2d(x, L, P, Q)

Computes autocorrelation wavelet coeficients for 2D signals.

# Arguments
- `x`: 2D signals (n ,m) matrix
- `L`: Decomposition level
- `P`: Low AC shell filter
- `Q`: High AC shell filter
"""
function ac2d(x,L,P,Q)
   onedim = ac1d_col(x,L,P,Q)
   twodim = [ac1d_row(i,L,P,Q) for i in onedim]
   return twodim
end

end # module
