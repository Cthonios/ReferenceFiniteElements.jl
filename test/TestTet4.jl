@testset ExtendedTestSet "Tet4 implementation" begin
  for int_type in [Int32, Int64]
    for float_type in [Float32, Float64]
      e = ReferenceFE(Tet4(), 1, int_type, float_type) # q_degree doesn't matter for this test
      v_nodes = vertex_nodes(e)
      @test e.nodal_coordinates[:, v_nodes[1]] ≈ [0., 0., 0.]
      @test e.nodal_coordinates[:, v_nodes[2]] ≈ [1., 0., 0.]
      @test e.nodal_coordinates[:, v_nodes[3]] ≈ [0., 1., 0.]
      @test e.nodal_coordinates[:, v_nodes[4]] ≈ [0., 0., 1.]
    end
  end
end  

common_test_sets(Tet4(), [1, 2], [Int32, Int64], [Float32, Float64])
