@testset "Wiggle" begin
  x = randn(8)
  y = acwt(x, wavelet(WT.db4))         # Normal ACDWT
  @test typeof(wiggle(y)) == Plots.Plot{Plots.GRBackend}
  @test typeof(wiggle!(y)) == Plots.Plot{Plots.GRBackend}
end