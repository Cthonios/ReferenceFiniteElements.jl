@testset ExtendedTestSet "ElementStencils.jl - Tri6 implementation" begin
  @testset "Test Tri6 element interpolant points in element" begin
    for q_degree in [1, 2, 3, 4, 5, 6]
      e = ElementStencil(Tri6(), q_degree)

      @test e.coordinates[:, e.vertex_nodes[1]] ≈ [0., 0.]
      @test e.coordinates[:, e.vertex_nodes[2]] ≈ [1., 0.]
      @test e.coordinates[:, e.vertex_nodes[3]] ≈ [0., 1.]
    end
  end  

  # @testset "Test Tri6 element face nodes match 1D lobatto nodes" begin
  #   for q_degree in [1, 2, 3, 4, 5, 6]
  #     e_1d = ElementStencil(Edge(), q_degree)
  #     e_2d = ElementStencil(Tri6(), q_degree)
  #     for face_node_ids in eachcol(e_2d.face_nodes)
  #       X_f = e_2d.coordinates[:, face_node_ids]
  #       temp = (1. .- e_1d.coordinates) * X_f[:, 1]' + e_1d.coordinates * X_f[:, end]'
  #       @test X_f ≈ temp'
  #     end 
  #   end
  # end

  @testset "Test Tri6 partition of unity property" begin
    for q_degree in [1, 2, 3, 4, 5, 6]
      partition_of_unity_tests(Tri6(), q_degree)
    end
  end

  @testset "Test Tri6 Kronecker delta property" begin
    for q_degree in [1, 2, 3, 4, 5, 6]
      kronecker_delta_property(Tri6(), q_degree)
    end
  end

  @testset "Test Tri6 interpolation" begin
    x = generate_random_points_in_triangle(1)
    q_degree = 1
    poly_coeffs = UpperTriangular(ones(Float64, q_degree + 1, q_degree + 1))
    poly_coeffs = map(x -> reverse(x), eachrow(poly_coeffs))
    poly_coeffs = mapreduce(permutedims, vcat, poly_coeffs)
    expected = polyval2d(x[1], x[2], poly_coeffs)

    e = ElementStencil(Tri6(), q_degree)
    Ns = ReferenceFiniteElements.shape_function_values(Tri6(), x)
    fn = polyval2d.(e.coordinates[1, :], e.coordinates[2, :], (poly_coeffs,))

    finterpolated = dot(Ns, fn)

    @test finterpolated ≈ expected
  end

  @testset "Test Tri6 grad interpolation" begin
    x = generate_random_points_in_triangle(1)
    q_degree = 1
    poly_coeffs = UpperTriangular(ones(Float64, q_degree + 1, q_degree + 1))
    poly_coeffs = map(x -> reverse(x), eachrow(poly_coeffs))
    poly_coeffs = mapreduce(permutedims, vcat, poly_coeffs)

    expected_dx = dpolyval2d(x[1], x[2], poly_coeffs, 1)
    expected_dy = dpolyval2d(x[1], x[2], poly_coeffs, 2)

    e = ElementStencil(Tri6(), q_degree)
    ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(Tri6(), x)
    fn = polyval2d.(e.coordinates[1, :], e.coordinates[2, :], (poly_coeffs,))

    temp_x = dot(∇N_ξ[:, 1], fn)
    temp_y = dot(∇N_ξ[:, 2], fn)

    @test temp_x ≈ expected_dx
    @test temp_y ≈ expected_dy
  end
end
