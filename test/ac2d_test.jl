using Plots, Wavelets, LinearAlgebra, FileIO, Images, BenchmarkTools

include("../src/AutocorrelationShell.jl")
using Main.AutocorrelationShell

img = load("./test/pictures/lenna.jpg")
img = Float64.(Gray.(img))

H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

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

function ac2d_matrix(x, L, P, Q)
    num_row, num_col = size(x)
    J = trunc(Integer, log2(num_col))
    D = J - L + 1
    accoef_matrix_3d = ac1d_col(img, L, P, Q)
    accoef_matrix_4d = Array{Float64, 4}(undef, D, num_row, D, num_col)
    for i in 1:D
        accoef_matrix_4d[i,:,:,:] = ac1d_col(accoef_matrix_3d[:,i,:],L,P,Q)
    end
    accoef_matrix_4d = permutedims(accoef_matrix_4d, [4,2,1,3])
    return accoef_matrix_4d
end

ac2d_matrix(img, L, P, Q)

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
        b = @benchmark ac2d_matrix(img, $L, $P, $Q)
        mean_runtime[idx] = mean(b).time
    end
    println(string("The average runtime was: "), mean(mean_runtime))
end

fwt_ac_benchmark("./test/pictures")
