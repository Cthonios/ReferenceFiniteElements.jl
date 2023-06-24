@testset ExtendedTestSet "ElementStencils.jl - Hex8 implementation" begin
  @testset "Test Hex8 element interpolant points in element" begin
    for q_degree in [1, 2, 3, 4, 5, 6]
      e = ElementStencil(Hex8(), q_degree)

      @test e.coordinates[:, e.vertex_nodes[1]] ≈ [-1.0, -1.0, -1.0]
      @test e.coordinates[:, e.vertex_nodes[2]] ≈ [1.0, -1.0, -1.0]
      @test e.coordinates[:, e.vertex_nodes[3]] ≈ [1.0, 1.0, -1.0]
      @test e.coordinates[:, e.vertex_nodes[4]] ≈ [-1.0, 1.0, -1.0]
      @test e.coordinates[:, e.vertex_nodes[5]] ≈ [-1.0, -1.0, 1.0]
      @test e.coordinates[:, e.vertex_nodes[6]] ≈ [1.0, -1.0, 1.0]
      @test e.coordinates[:, e.vertex_nodes[7]] ≈ [1.0, 1.0, 1.0]
      @test e.coordinates[:, e.vertex_nodes[8]] ≈ [-1.0, 1.0, 1.0]
    end
  end  

  @testset "Test Quad4 element shape function partition of unity" begin
    for q_degree in [1, 2, 3, 4, 5, 6]
      partition_of_unity_tests(Hex8(), q_degree)
    end 
  end

  @testset "Test Quad4 element Kronecker delta property" begin
    for q_degree in [1, 2, 3, 4, 5, 6]
      kronecker_delta_property(Hex8(), q_degree)
    end
  end
end
