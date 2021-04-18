@testset "ACWD" begin
  x₁ = randn(8)
  y₁ = acwt(x₁, wavelet(WT.db4))
  @test iacwt(y₁) ≈ x₁

  x₂ = randn(8,8)
  y₂ = acwt(x₂, wavelet(WT.db4))
  @test iacwt(y₂) ≈ x₂
end

@testset "ACWPT" begin
  x₁ = randn(8)
  y₁ = acwt(x₁, wavelet(WT.db4))
  y₂ = acwpt(x₁, wavelet(WT.db4))
  
  @test y₂[:,8] ≈ y₁[:,1]
  @test iacwpt(y₂, maketree(x₁, :full)) ≈ x₁
end
