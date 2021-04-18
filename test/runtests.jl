using 
    Test,
    AutocorrelationShell, 
    Wavelets, 
    LinearAlgebra, 
    Plots

@testset "Transforms" begin include("testsets/transforms.jl") end
@testset "Visualizations" begin include("testsets/visualizations.jl") end

