@testset ExtendedTestSet "Tet10 implementation" begin
  for int_type in [Int32, Int64]
    for float_type in [Float32, Float64]
      for array_type in [SArray]
        e = ReferenceFE(Tet10(Val(1)); int_type=int_type, float_type=float_type, array_type=array_type) # q_degree doesn't matter for this test
        v_nodes = vertex_nodes(e)
        @test e.nodal_coordinates[:, v_nodes[1]]  ≈ [0., 0., 0.]
        @test e.nodal_coordinates[:, v_nodes[2]]  ≈ [1., 0., 0.]
        @test e.nodal_coordinates[:, v_nodes[3]]  ≈ [0., 1., 0.]
        @test e.nodal_coordinates[:, v_nodes[4]]  ≈ [0., 0., 1.]
        #
        @test e.nodal_coordinates[:, v_nodes[5]]  ≈ [0.5, 0.0, 0.0]
        @test e.nodal_coordinates[:, v_nodes[6]]  ≈ [0.5, 0.5, 0.0]
        @test e.nodal_coordinates[:, v_nodes[7]]  ≈ [0.0, 0.5, 0.0]
        @test e.nodal_coordinates[:, v_nodes[8]]  ≈ [0.0, 0.0, 0.5]
        @test e.nodal_coordinates[:, v_nodes[9]]  ≈ [0.5, 0.0, 0.5]
        @test e.nodal_coordinates[:, v_nodes[10]] ≈ [0.0, 0.5, 0.5]
      end
    end
  end
end  

# erroring on gradient partition of unity check
common_test_sets(Tet10, [1, 2], [Int32, Int64], [Float32, Float64], [SArray], [StructArray])

if CUDA.has_cuda()
  common_test_sets(Tet10, [1, 2], [Int32, Int64], [Float32, Float64], [SArray], [StructArray]; cuda=true)
end
