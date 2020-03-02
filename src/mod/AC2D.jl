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
