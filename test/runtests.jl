using Test, AutocorrelationShell, Wavelets, LinearAlgebra, AbstractTrees, Random

@testset "Autocorrelation Shell" begin
    Q = qfilter(wavelet(WT.db2));
    P = pfilter(wavelet(WT.db2));
    x = zeros(256); x[128] = 1;

    # Autocorrelation Wavelet Transform 1D
    @testset "ACW1D" begin
        @test begin
            decomp = acwt(x, L=0, P=P, Q=Q)
            norm(x - iacwt(decomp)) < 1e-15
        end
    end

    # Autocorrelation Wavelet Transform 2D
    @testset "ACW2D" begin
        @test begin
            X = randn(512, 512);
            decomp = acwt2D(X; L_row=4, L_col=4, P=P, Q=Q)
            norm(X - iacwt2D(decomp)) < 1e-12
        end
    end

    # Autocorrelation Wavelet Packet Transform
    @testset "ACWPT" begin
        @test begin
            X = acwt(x, L=2, P=P, Q=Q)[:,4];
            decomp = acwt(X, L=0, P=P, Q=Q);
            tree = acwpt(X, P, Q);
            norm(collect(PostOrderDFS(tree))[1].data - decomp[:,1]) == 0
        end

        @test begin
            X = acwt(x, L=2, P=P, Q=Q)[:,4];
            tree = acwpt(X, P, Q);
            norm(X - iacwpt(tree)) < 1e-14
        end

        @test begin
            X = acwt(x, L=2, P=P, Q=Q)[:,4];
            tree = acwpt(X, P, Q);
            best_tree = acwptPostOrderBestBasis(tree)
            norm(X - iacwpt(best_tree)) < 1e-15
        end

        @test begin
            X = acwt(x, L=2, P=P, Q=Q)[:,4];
            tree = acwpt(X, P, Q);
            best_tree = acwptPreOrderBestBasis(tree)
            norm(X - iacwpt(best_tree)) < 1e-15
        end
    end

    @testset "ACUtil" begin
        @test snr([0, 100], [1, 100]) == 40

        @test acthreshold([0.1, 0.2, 0.2, 0.3], "hard", 0.1) == [0, 0.2, 0.2, 0.3]
        @test acthreshold([0.05, 0.2, 0.4, 0.5], "soft", 0.1) ≈ [0, 0.1, 0.3, 0.4]
    end
end
