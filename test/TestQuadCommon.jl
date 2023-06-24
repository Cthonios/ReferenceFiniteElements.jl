@testset ExtendedTestSet "Quadratures.jl - QuadCommon implementation" begin
  @testset "Test quad weight positivity" begin
    for q_degree in [1, 2, 3, 4, 5, 6]
      q_rule = Quadrature(Quad4(), q_degree)
      for w in q_rule.ws
        @test w > 0.
      end
    end
  end

  @testset "Test sum of quadrature points" begin
    for q_degree in [1, 2, 3, 4, 5, 6]
      q_rule = Quadrature(Quad4(), q_degree)
      @test q_rule.ws |> sum ≈ 4.0
    end
  end

  @testset "Test quad point in domain" begin
    for q_degree in [1, 2, 3, 4, 5, 6]
      q_rule = Quadrature(Quad4(), q_degree)
      for ξ in eachcol(q_rule.ξs)
        @test is_inside_quad(ξ)
      end
    end
  end
end

