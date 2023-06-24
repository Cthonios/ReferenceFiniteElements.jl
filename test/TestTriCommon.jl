@testset ExtendedTestSet "Quadratures.jl - TriCommon implementation" begin
  @testset "Test quad weight positivity" begin
    for degree in [1, 2, 3, 4, 5, 6]
      q_rule = Quadrature(Tri3(), degree)
      for w in q_rule.ws
        @test w > 0.
      end
    end
  end

  @testset "Test quad point in domain" begin
    for degree in [1, 2, 3, 4, 5, 6]
      q_rule = Quadrature(Tri3(), degree)
      for ξ in eachcol(q_rule.ξs)
        @test is_inside_triangle(ξ)
      end
    end
  end

  @testset "Test triangle quadrature exactness" begin
    for degree in [1, 2, 3, 4, 5, 6]
      q_rule = Quadrature(Tri3(), degree)
      for i in 1:degree
        for j in 1:degree - i
          quad_answer = map((ξ, w) -> w * ξ[1]^i * ξ[2]^j, eachcol(q_rule.ξs), q_rule.ws) |> sum
          exact = integrate_2d_monomial_on_triangle(i, j)
          @test quad_answer ≈ exact
        end
      end
    end
  end
end
