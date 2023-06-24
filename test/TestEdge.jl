end_points = [0., 1.]
max_degree_1d = 25
max_degree = 6

@testset ExtendedTestSet "ElementStencils.jl - Edge implementation" begin
  @testset "Test 1D interpolant points in element" begin
    for degree in 1:max_degree
      e = ElementStencil(Edge(), degree)
      @test e.degree == degree
      @test all(X -> X >= 0. && X <= 1., e.coordinates)
    end
  end

  @testset "Test 1D element topological nodesets" begin
    for degree in 1:max_degree
      e = ElementStencil(Edge(), degree)
      @test e.coordinates[e.vertex_nodes[1]] ≈ 0.
      @test e.coordinates[e.vertex_nodes[2]] ≈ 1.

      if length(e.interior_nodes) > 0
        @test all(x -> x > 0. && x < 1., e.coordinates[e.interior_nodes])
      end
    end
  end
end

@testset ExtendedTestSet "Quadratures.jl - Edge implementation" begin
  @testset "Test 1d quad weight positivity" begin
    for degree in 1:max_degree_1d
      q_rule = Quadrature(Edge(), degree)
      @test all(w -> w > 0., q_rule.ws)
    end
  end

  @testset "Test 1d quadrature points in domain" begin
    for degree in 1:max_degree_1d
      q_rule = Quadrature(Edge(), degree)
      @test all(ξ -> is_inside_unit_interval(ξ), q_rule.ξs)
    end
  end

  # TODO fix to be like optimisim
  # there's something weird in the 1d shape functions there
  @testset "Test 1d quadrature exactness" begin
    for degree in 1:max_degree_1d
      q_rule = Quadrature(Edge(), degree)
      for i in 1:degree
        x_vals = map_affine_1d.(q_rule.ξs, (end_points,))
        jac = map_1d_jac(end_points)
        monomial = x_vals.^i
        quad_answer = sum(dot(monomial, q_rule.ws)) * jac
        exact_answer = integrate_1d_monomial_on_line(i)
        @test quad_answer ≈ exact_answer
      end
    end
  end
end
