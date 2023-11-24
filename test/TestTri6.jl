@testset ExtendedTestSet "Tri6 implementation" begin
  for int_type in [Int32, Int64]
    for float_type in [Float32, Float64]
      e = ReferenceFE(Tri6(1), int_type, float_type)
      v_nodes = vertex_nodes(e)
      @test e.nodal_coordinates[:, v_nodes[1]] ≈ [0., 0.]
      @test e.nodal_coordinates[:, v_nodes[2]] ≈ [1., 0.]
      @test e.nodal_coordinates[:, v_nodes[3]] ≈ [0., 1.]
      @test e.nodal_coordinates[:, v_nodes[4]] ≈ [0.5, 0.0]
      @test e.nodal_coordinates[:, v_nodes[5]] ≈ [0.5, 0.5]
      @test e.nodal_coordinates[:, v_nodes[6]] ≈ [0.0, 0.5]
    end
  end
end

# @testset "Test Tri6 element face nodes match 1D lobatto nodes" begin
#   for q_degree in [1, 2, 3, 4, 5, 6]
#     e_1d = ReferenceFEStencil(Edge(), q_degree)
#     e_2d = ReferenceFEStencil(Tri6(), q_degree)
#     for face_node_ids in eachcol(e_2d.face_nodes)
#       X_f = e_2d.coordinates[:, face_node_ids]
#       temp = (1. .- e_1d.coordinates) * X_f[:, 1]' + e_1d.coordinates * X_f[:, end]'
#       @test X_f ≈ temp'
#     end 
#   end
# end


@testset ExtendedTestSet "Test Tri6 interpolation" begin
  for int_type in [Int32, Int64]
    for float_type in [Float32, Float64]
      x = generate_random_points_in_triangle(1)
      x = reinterpret(SVector{2, float_type}, vec(x))
      x = x[1]
      q_degree = 1
      poly_coeffs = UpperTriangular(ones(float_type, q_degree + 1, q_degree + 1))
      poly_coeffs = map(x -> reverse(x), eachrow(poly_coeffs))
      poly_coeffs = mapreduce(permutedims, vcat, poly_coeffs)
      expected = polyval2d(x[1], x[2], poly_coeffs)
    
      e = ReferenceFE(Tri6(q_degree), int_type, float_type)
      # Ns = ReferenceFiniteElements.shape_function_values(Tri6(q_degree), x)
      # TODO currently defaulting to SVector
      Ns = ReferenceFiniteElements.shape_function_values(Tri6(q_degree), SVector, x)
      fn = polyval2d.(e.nodal_coordinates[1, :], e.nodal_coordinates[2, :], (poly_coeffs,))
    
      finterpolated = dot(Ns, fn)
    
      if float_type == Float32
        @test_skip finterpolated ≈ expected atol=1e-7 rtol=1e-7
      elseif float_type == Float64
        @test finterpolated ≈ expected
      end
    end
  end
end

@testset ExtendedTestSet "Test Tri6 grad interpolation" begin
  for int_type in [Int32, Int64]
    for float_type in [Float32, Float64]
      x = generate_random_points_in_triangle(1)
      x = reinterpret(SVector{2, float_type}, vec(x))
      x = x[1]
      q_degree = 1
      poly_coeffs = UpperTriangular(ones(Float64, q_degree + 1, q_degree + 1))
      poly_coeffs = map(x -> reverse(x), eachrow(poly_coeffs))
      poly_coeffs = mapreduce(permutedims, vcat, poly_coeffs)

      expected_dx = dpolyval2d(x[1], x[2], poly_coeffs, 1)
      expected_dy = dpolyval2d(x[1], x[2], poly_coeffs, 2)

      e = ReferenceFE(Tri6(q_degree), int_type, float_type)
      # ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(Tri6(q_degree), x)
      # TODO defaulting to SMatrix for now
      ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(Tri6(q_degree), SMatrix, x)
      fn = polyval2d.(e.nodal_coordinates[1, :], e.nodal_coordinates[2, :], (poly_coeffs,))

      temp_x = dot(∇N_ξ[:, 1], fn)
      temp_y = dot(∇N_ξ[:, 2], fn)

      if float_type == Float32
        @test_skip temp_x ≈ expected_dx atol=1e-7 rtol=1e-7
        @test_skip temp_y ≈ expected_dy atol=1e-7 rtol=1e-7
      elseif float_type == Float64
        @test temp_x ≈ expected_dx
        @test temp_y ≈ expected_dy
      end
    end
  end
end

common_test_sets(Tri6, [1, 2, 3, 4, 5, 6], [Int32, Int64], [Float32, Float64])
