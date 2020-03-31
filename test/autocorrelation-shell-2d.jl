# Add packages
using Pkg, Plots, Wavelets, LinearAlgebra, Statistics, FileIO, Images, BenchmarkTools

# AutocorrelationShell -------------------
# Include the AutocorrelationShell module and wiggle function
include("./wiggle.jl")
include("./src/AutocorrelationShell.jl")
include("../src/ACShell.jl")
using Main.ACShell

# Instanciate Wavelets
H = wavelet(WT.db2)
L = 4
Q = qfilter(H)
P = pfilter(H)

img = load("./pictures/lenna.jpg")
img = Float64.(Gray.(img))
heatmap(img)

# One dimensional functions ------------------
"""
    ac1d_col(x, L, P, Q)

Computes the column wise one dimension fwt_ac coeficients.

# Arguments
- `x`: Signal(Image)
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

# two dimensional functions ------------------
function ac2d_rowcol(x,L,P,Q)
    onedim = ac1d_row(x,L,P,Q)
    twodim = [ac1d_col(i,L,P,Q) for i in onedim]
    return twodim
end

function ac2d_colrow(x,L,P,Q)
    onedim = ac1d_col(x,L,P,Q)
    twodim = [ac1d_row(i,L,P,Q) for i in onedim]
    return twodim
end


fwt_img = fwt_ac_twodim(img)
size(fwt_img)


"""
    ac2d_heatmap(x)

Visualizes the result of fwt_ac_twodim.

# Arguments
- `ac2d`: output of ac2d function i.e. vector of coefficients
"""
function ac2d_heatmap(ac2d)
    x = [reduce(hcat, ac2d) for i in 1:size(ac2d)]
    x = reduce(vcat, x)
    x = abs.(x) # Ensure that there are no negative values
    x = log.(x)
    heatmap(x)
end

# Compute L2 Norm for rowcol and colrow implementation
# NOTE: the result of rowcol is the transpose of colrow
rowcol = ac2d_rowcol(img,L,P,Q)
colrow = ac2d_colrow(img,L,P,Q)

function ac_re(x,y)
    s=size(x)[1]
    re=Array{Float64}(undef,s,s)
    for i in 1:s
        for j in 1:s
            re[i,j] = norm(x[i][j][:] - y[j][i][:])/norm(x[i][j][:])
        end
    end
    return re
end

mean(ac_re(rowcol,colrow))

# Benchmarking ----------------
"""
    fwt_ac_benchmark(PATH, f)

Computes the mean of average runtime for fwt_ac_twodim on test images.

# Arguments
- `PATH`: direct or relative path to the directory containing test images.
"""
function fwt_ac_benchmark(PATH)
    fn_list = [fn for fn in readdir(PATH) if fn != ".DS_Store"]
    mean_runtime = Vector{Float64}(undef, length(fn_list))
    for (idx, fn) in enumerate(fn_list)
        img = load(string(PATH, "/", fn))
        println(string("Computing runtime for: ", fn))
        img = Float32.(Gray.(img))
        b = @benchmark ac2d(img, $L, $P, $Q)
        mean_runtime[idx] = mean(b).time
    end
    println(string("The average runtime was: "), mean(mean_runtime))
end

fwt_ac_benchmark("./test/pictures")
