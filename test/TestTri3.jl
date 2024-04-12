#   # @testset "Test Tri3 element face nodes match 1D lobatto nodes" begin
#   #   for q_degree in [1, 2, 3, 4, 5, 6]
#   #     e_1d = ReferenceFEStencil(Edge(), q_degree)
#   #     e_2d = ReferenceFEStencil(Tri3(), q_degree)
#   #     for face_node_ids in eachcol(e_2d.face_nodes)
#   #       X_f = e_2d.coordinates[:, face_node_ids]
#   #       temp = (1. .- e_1d.coordinates) * X_f[:, 1]' + e_1d.coordinates * X_f[:, end]'
#   #       @test X_f ≈ temp'
#   #     end 
#   #   end
#   # end

@testset ExtendedTestSet "Tri3 implementation" begin
  for int_type in [Int32, Int64]
    for float_type in [Float32, Float64]
      for array_type in [SArray]
        e = ReferenceFE(Tri3(Val(1)); int_type=int_type, float_type=float_type, array_type=array_type)
        v_nodes = vertex_nodes(e)
        @test e.nodal_coordinates[:, v_nodes[1]] ≈ [0., 0.]
        @test e.nodal_coordinates[:, v_nodes[2]] ≈ [1., 0.]
        @test e.nodal_coordinates[:, v_nodes[3]] ≈ [0., 1.]
      end
    end
  end
end

# @testset ExtendedTestSet "Tri3 element face nodes match 1D lobatto nodes" begin
#   e_1d = ReferenceFE(Edge(), 1)
#   e_2d = ReferenceFE(Tri3(), 1)
#   for face_node_ids in eachcol(e_2d.face_nodes)
#     X_f = e_2d.nodal_coordinates[:, face_node_ids]
#     temp = (1. .- e_1d.nodal_coordinates) * X_f[:, 1]' + e_1d.nodal_coordinates * X_f[:, end]'
#     @test X_f ≈ temp'
#   end 
# end

# @testset ExtendedTestSet "Test Tri3 partition of unity property" begin
#   for q_degree in [1, 2, 3, 4, 5, 6]
#     partition_of_unity_tests(Tri3(), q_degree)
#   end
# end

# @testset ExtendedTest"Test Tri3 Kronecker delta property" begin
#   for q_degree in [1, 2, 3, 4, 5, 6]
#     kronecker_delta_property(Tri3(), q_degree)
#   end
# end

@testset ExtendedTestSet "Test Tri3 interpolation" begin
  for int_type in [Int32, Int64]
    for float_type in [Float32, Float64]
      for array_type in [SArray]
        x = generate_random_points_in_triangle(1)
        x = reinterpret(SVector{2, float_type}, vec(x))
        x = x[1]
        q_degree = 1
        poly_coeffs = UpperTriangular(ones(float_type, q_degree + 1, q_degree + 1))
        poly_coeffs = map(x -> reverse(x), eachrow(poly_coeffs))
        poly_coeffs = mapreduce(permutedims, vcat, poly_coeffs)
        expected = polyval2d(x[1], x[2], poly_coeffs)
        e = ReferenceFE(Tri3(Val(1)); int_type=int_type, float_type=float_type, array_type=array_type)
        # Ns = ReferenceFiniteElements.shape_function_values(Tri3(1), x)
        # TODO currently defaulting to SVector in these tests
        Ns = ReferenceFiniteElements.shape_function_values(Tri3(Val(1)), SVector, x)
        fn = polyval2d.(e.nodal_coordinates[1, :], e.nodal_coordinates[2, :], (poly_coeffs,))
        finterpolated = dot(Ns, fn)
        if float_type == Float32
          @test finterpolated ≈ expected atol=1e-7 rtol=1e-7
        elseif float_type == Float64
          @test finterpolated ≈ expected
        end
      end
    end
  end
end

@testset ExtendedTestSet "Test Tri3 grad interpolation" begin
  for int_type in [Int32, Int64]
    for float_type in [Float32, Float64]
      for array_type in [SArray]
        x = generate_random_points_in_triangle(1)
        x = reinterpret(SVector{2, float_type}, vec(x))
        x = x[1]
        q_degree = 1
        poly_coeffs = UpperTriangular(ones(float_type, q_degree + 1, q_degree + 1))
        poly_coeffs = map(x -> reverse(x), eachrow(poly_coeffs))
        poly_coeffs = mapreduce(permutedims, vcat, poly_coeffs)

        expected_dx = dpolyval2d(x[1], x[2], poly_coeffs, 1)
        expected_dy = dpolyval2d(x[1], x[2], poly_coeffs, 2)
        e = ReferenceFE(Tri3(Val(q_degree)); int_type=int_type, float_type=float_type, array_type=array_type)
        # ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(Tri3(q_degree), x)
        # TODO defaulting to SMatrix for now
        ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(Tri3(Val(q_degree)), SMatrix, x)
        fn = polyval2d.(e.nodal_coordinates[1, :], e.nodal_coordinates[2, :], (poly_coeffs,))

        temp_x = dot(∇N_ξ[:, 1], fn)
        temp_y = dot(∇N_ξ[:, 2], fn)

        if float_type == Float32
          @test temp_x ≈ expected_dx atol=1e-7 rtol=1e-7
          @test temp_y ≈ expected_dy atol=1e-7 rtol=1e-7
        elseif float_type == Float64
          @test temp_x ≈ expected_dx
          @test temp_y ≈ expected_dy
        end
      end
    end
  end
end

common_test_sets(Tri3, [1, 2], [Int32, Int64], [Float32, Float64], [SArray], [StructArray])

if CUDA.has_cuda()
  common_test_sets(Tri3, [1, 2], [Int32, Int64], [Float32, Float64], [SArray], [StructArray]; cuda=true)
end
