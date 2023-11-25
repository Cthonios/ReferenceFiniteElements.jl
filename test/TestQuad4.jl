@testset ExtendedTestSet "Quad4 implementation" begin
  for int_type in [Int32, Int64]
    for float_type in [Float32, Float64]
      for array_type in [SArray, MArray]
        e = ReferenceFE(Quad4(1); int_type=int_type, float_type=float_type, array_type=array_type) # q_degree doesn't matter for this test
        v_nodes = vertex_nodes(e)
        coords = e.nodal_coordinates
        @test e.nodal_coordinates[:, v_nodes[1]] ≈ [-1.0, -1.0]
        @test e.nodal_coordinates[:, v_nodes[2]] ≈ [1.0, -1.0]
        @test e.nodal_coordinates[:, v_nodes[3]] ≈ [1.0, 1.0]
        @test e.nodal_coordinates[:, v_nodes[4]] ≈ [-1.0, 1.0]
      end
    end
  end
end  

common_test_sets(Quad4, [1, 2, 3, 4, 5, 6], [Int32, Int64], [Float32, Float64], [SArray, MArray])
