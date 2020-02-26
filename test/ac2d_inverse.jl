# --- Implement inverse ac2d ---

# Import packages
using Pkg, Plots, Wavelets, LinearAlgebra, FileIO, Images

include("../src/AutocorrelationShell.jl")
using Main.AutocorrelationShell

# Load test image
img = load("./test/pictures/lenna.jpg")
img = img = Float64.(Gray.(img))

# Instanciate Wavelets
H = wavelet(WT.db2)
L = 4
Q = qfilter(H)
P = pfilter(H)

# Perform decomposition
decomp = ac2d(img, L, P, Q)

# Implement inverse ac2d
function inv_ac2d(ac2d)
    """
        inv_ac2d(ac2d)

        Recreates the original 2D signal from the coeficients produced by ac2d function.

        # Arguments
        `ac2d`: The output of ac2d i.e. vector autocorrelation wavelet transform coefficients.
    """
    n = size(ac2d)[1]
    y = deepcopy(ac2d[n][n]) # Deepcopy finest scale
    for i in 1:n-1
        yi = deepcopy(ac2d[i][n])
        for j in 1:n-1
            yi = (yi + ac2d[i][j]) / sqrt(2)
        end
        y = (y + yi) / sqrt(2)
    end
    return y
end

reconst = inv_ac2d(decomp)

# Compute difference between original image
norm(reconst-img) / norm(img)

# Compare visually
heatmap(img)
heatmap(reconst)
