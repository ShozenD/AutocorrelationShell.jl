using Test, AutocorrelationShell, Wavelets, LinearAlgebra, Plots

x = zeros(256); x[128] = 1;
X = acwt(x, wavelet(WT.db4))[:,4];
y₁ = acwpt(X, wavelet(WT.db4)); # Array method
bb = bestbasistree(y₁, BB(stationary = true))

@testset "Autocorrelation Shell" begin
    x = zeros(256); x[128] = 1;

    # Autocorrelation Wavelet Transform 1D
    @testset "ACTransforms" begin

        @testset "1D transform" begin
            y = acwt(x, wavelet(WT.db4))
            @test norm(x - iacwt(y)) < 1e-15
        end

        @testset "2D transform" begin
            X = randn(256, 256);
            decomp = acwt(X, wavelet(WT.db4))
            @test norm(X - iacwt(decomp)) < 1e-12
        end
    end

    # Autocorrelation Wavelet Packet Transform
    @testset "ACWPT" begin
        X = acwt(x, wavelet(WT.db4))[:,4];

        y₁ = acwpt(X, wavelet(WT.db4)); # Array method
        y = acwt(X, wavelet(WT.db4)); # Normal ACDWT
        
        @test norm(y₁[:,256] - y[:,1]) < 1e-15

        bb = bestbasistree(y₁, BB(stationary = true))
        @test norm(X - iacwpt(y₁,bb)) < 1e-15

        # p = selectednodes_plot(besttree) 
        # @test typeof(p) == Plots.Plot{Plots.GRBackend}

        # p = selectednodes_plot(bb)
        # @test typeof(p) == Plots.Plot{Plots.GRBackend}
    end

    @testset "ACPlots" begin
        y = acwt(x, wavelet(WT.db4)); # Normal ACDWT
        @test typeof(wiggle(y)) == Plots.Plot{Plots.GRBackend}
        @test typeof(wiggle!(y)) == Plots.Plot{Plots.GRBackend}
    end

    @testset "ACUtil" begin
        X = randn(128,128)
        y = acwt(X, wavelet(WT.db4)); 
        @test typeof(acwt_heatmap(y[1,1,:,:])) == Plots.Plot{Plots.GRBackend}
    end
end