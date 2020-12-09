using
    Test,
    AutocorrelationShell,
    Wavelets,
    LinearAlgebra,
    AbstractTrees,
    Random,
    Plots,
    CSV,
    DataFrames

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
        X = acwt(x, L=2, P=P, Q=Q)[:,4];
        tree = acwpt(X, P, Q);

        @test begin
            decomp = acwt(X, L=0, P=P, Q=Q);
            norm(collect(PostOrderDFS(tree))[1].data - decomp[:,1]) == 0
        end

        @test begin
            norm(X - iacwpt(tree)) < 1e-14
        end

        @test begin
            best_tree = acwpt_postorder_bb(tree)
            norm(X - iacwpt(best_tree)) < 1e-15
        end

        @test begin
            best_tree = acwpt_preorder_bb(tree)
            norm(X - iacwpt(best_tree)) < 1e-15
        end

        @test begin
            best_tree = acwpt_postorder_bb(tree, et=ShannonEntropy())
            norm(X - iacwpt(best_tree)) < 1e-15
        end

        @test begin
            best_tree = acwpt_postorder_bb(tree, et=LogEnergyEntropy())
            norm(X - iacwpt(best_tree)) < 1e-15
        end

        @test begin
            best_tree = acwpt_postorder_bb(tree)
            p = selectednodes_plot(best_tree)
            typeof(p) == Plots.Plot{Plots.GRBackend}
        end
    end

    @testset "ACUtil" begin
        @test snr([0, 100], [1, 100]) == 40

        @test acthreshold([0.1, 0.2, 0.2, 0.3], "hard", 0.1) == [0, 0.2, 0.2, 0.3]
        @test acthreshold([0.05, 0.2, 0.4, 0.5], "soft", 0.1) ≈ [0, 0.1, 0.3, 0.4]

        # For the sake of code coverage
        @test begin
            p = wiggle(acwt(x, L=1, P=P, Q=Q))
            typeof(p) == Plots.Plot{Plots.GRBackend}

            p = wiggle!(acwt(x, L=1, P=P, Q=Q))
            typeof(p) == Plots.Plot{Plots.GRBackend}
        end

        # Test denoising
        @test begin
            test_data = CSV.File("data/wavelet_test_256.csv") |> DataFrame
            rng = MersenneTwister(123);
            y = test_data.doppler;
            y_noisy = make_noisy(y, rng, 0.07);

            acwpt_decomp = acwpt(y_noisy, P, Q);
            bb = acwpt_postorder_bb(acwpt_decomp);

            acwt_decomp = acwt(y_noisy; L=1, P=P, Q=Q);

            t1, s1 = acwt_snr(y, y_noisy, "soft", 0.01);
            thresh1 = acthreshold(acwt_decomp, "soft", t1[findmax(s1)[2]]);
            acwt_reconst = iacwt(thresh1);

            t2, s2 = acwpt_snr(y, y_noisy, "soft", 0.01);
            thresh2 = threshold_bestbasis(bb, "soft", t2[findmax(s2)[2]]);
            acwpt_reconst = iacwpt(thresh2);

            (typeof(acwt_reconst) == Vector{Float64}) & (typeof(acwpt_reconst) == Vector{Float64})
        end
    end
end
