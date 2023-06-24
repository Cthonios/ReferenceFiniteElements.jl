@testset ExtendedTestSet "Quadratures.jl - TetCommon implementation" begin
  @testset "Test tet point in domain" begin
    for degree in [1, 2]
      q_rule = Quadrature(Tet4(), degree)
      for ξ in eachcol(q_rule.ξs)
        @test is_inside_tet(ξ)
      end
    end
  end

  @testset "Test tet quadrature weight sums" begin
    for degree in [1, 2]
      q_rule = Quadrature(Tet4(), degree)
      @test q_rule.ws |> sum ≈ 1. / 6.
    end
  end
end
