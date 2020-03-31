# --- Implement inverse ac2d ---
# Import packages
using Pkg, Plots, Wavelets, LinearAlgebra, FileIO, Images

include("../src/AutocorrelationShell.jl")
using Main.AutocorrelationShell

# Load test image
img = load("./test/pictures/lenna.jpg")
img = Float64.(Gray.(img))

# Instanciate Wavelets
H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

# Perform decomposition
decomp = ac2d(img, L, P, Q)

function iac2d_matrix(x)
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

reconst = iac2d_matrix(decomp)

# Compute difference between original image
norm(reconst[:]-img[:]) / norm(img[:])

# Compare visually
heatmap(img)
heatmap(reconst)
