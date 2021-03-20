using Test, AutocorrelationShell, Wavelets, LinearAlgebra, Plots

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
        y₁ = acwpt(X, wavelet(WT.db4));
        y₂ = acwt(X, wavelet(WT.db4));
        @test norm(y₁[:,256] - y₂[:,1]) < 1e-15

        bb = bestbasistree(y₁, NormEntropy())
        @test norm(X - iacwpt(y₁,bb)) < 1e-15

        bb = bestbasistree(y₁, ShannonEntropy())
        @test norm(X - iacwpt(y₁,bb)) < 1e-15

        bb = bestbasistree(y₁, LogEnergyEntropy())
        @test norm(X - iacwpt(y₁,bb)) < 1e-15

        bb = bestbasistree(y₁, NormEntropy())
        p = selectednodes_plot(bb)
        @test typeof(p) == Plots.Plot{Plots.GRBackend}
    end
end
