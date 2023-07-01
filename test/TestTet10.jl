@testset ExtendedTestSet "ElementStencils.jl - Tet10 implementation" begin
  @testset "Test Tet10 element interpolant points in element" begin
    for q_degree in [1, 2]
      e = ElementStencil(Tet10(), q_degree)

      @test e.coordinates[:, e.vertex_nodes[1]] ≈ [0., 0., 0.]
      @test e.coordinates[:, e.vertex_nodes[2]] ≈ [1., 0., 0.]
      @test e.coordinates[:, e.vertex_nodes[3]] ≈ [0., 1., 0.]
      @test e.coordinates[:, e.vertex_nodes[4]] ≈ [0., 0., 1.]
    end
  end 

  @testset "Test Tet10 partition of unity property" begin
    for q_degree in [1, 2]
      partition_of_unity_tests(Tet4(), q_degree)
    end
  end

  @testset "Test Tet10 Kronecker delta property" begin
    for q_degree in [1, 2]
      kronecker_delta_property(Tet4(), q_degree)
    end
  end
end