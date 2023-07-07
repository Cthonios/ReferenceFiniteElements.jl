# max_degree = 6

# @testset ExtendedTestSet "ReferenceFEStencils.jl - SimplexTri implementation" begin
#   @testset "Test SimplexTri element interpolant points in element" begin
#     for degree in 1:max_degree
#       e = ReferenceFEStencil(SimplexTri(), degree)

#       @test e.coordinates[:, e.vertex_nodes[1]] ≈ [1., 0.]
#       @test e.coordinates[:, e.vertex_nodes[2]] ≈ [0., 1.]
#       @test e.coordinates[:, e.vertex_nodes[3]] ≈ [0., 0.]

#       if length(e.interior_nodes) > 0
#         @test all(p -> p[1]             > -eps(Float64), eachcol(e.coordinates[:, e.interior_nodes]))
#         @test all(p -> p[2] + p[1] - 1. <  eps(Float64), eachcol(e.coordinates[:, e.interior_nodes]))
#       end
#     end
#   end

#   @testset "Test SimplexTri element face nodes match 1D lobatto nodes" begin
#     for degree in 1:max_degree
#       e_1d = ReferenceFEStencil(Edge(), degree)
#       e_2d = ReferenceFEStencil(SimplexTri(), degree)
#       for face_node_ids in eachcol(e_2d.face_nodes)
#         X_f = e_2d.coordinates[:, face_node_ids]
#         X_f_1d = (1. .- e_1d.coordinates) * X_f[:, 1]' + e_1d.coordinates * X_f[:, end]'
#         @test X_f ≈ X_f_1d'
#       end
#     end
#   end

#   @testset "Test SimplexTri element shape function values partition of unity" begin
#     for degree in 1:max_degree
#       partition_of_unity_shape_function_values_int_test(SimplexTri(), degree)
#     end
#   end  

#   @testset "Test SimplexTri element shape function gradients partition of unity" begin
#     for degree in 1:max_degree
#       partition_of_unity_shape_function_gradients_int_test(SimplexTri(), degree)
#     end
#   end

#   @testset "Test SimplexTri element shape kronecker delta property" begin
#     for degree in 1:max_degree
#       q_rule = Quadrature(SimplexTri(), degree)
#       e = ReferenceFEStencil(SimplexTri(), degree)
#       n_nodes = (degree + 1) * (degree + 2) ÷ 2

#       A  = zeros(Float64, n_nodes, size(e.coordinates, 2))
#       Ax = zeros(Float64, n_nodes, size(e.coordinates, 2))
#       Ay = zeros(Float64, n_nodes, size(e.coordinates, 2))

#       nf  = zeros(Float64, n_nodes, size(e.coordinates, 2))
#       nfx = zeros(Float64, n_nodes, size(e.coordinates, 2))
#       nfy = zeros(Float64, n_nodes, size(e.coordinates, 2))

#       ReferenceFEs.vander2d!(A, Ax, Ay, e.coordinates, degree)
#       ReferenceFEs.vander2d!(nf, nfx, nfy, e.coordinates, degree)

#       Ns = A \ nf
#       @test Ns ≈ I
#     end
#   end

#   @testset "Test SimplexTri interpolation" begin
#     x = generate_random_points_in_triangle(1)
#     for degree in 1:max_degree
#       e = ReferenceFEStencil(SimplexTri(), degree)
#       n_nodes = (degree + 1) * (degree + 2) ÷ 2

#       # poly_coeffs = UpperTriangular(ones(Float64, degree + 1, degree + 1))
#       poly_coeffs = UpperTriangular(ones(Float64, degree + 1, degree + 1))
#       poly_coeffs = map(x -> reverse(x), eachrow(poly_coeffs))
#       poly_coeffs = mapreduce(permutedims, vcat, poly_coeffs)
#       expected = polyval2d(x[1], x[2], poly_coeffs)

#       # get shape functions
#       A  = zeros(Float64, n_nodes, size(e.coordinates, 2))
#       Ax = zeros(Float64, n_nodes, size(e.coordinates, 2))
#       Ay = zeros(Float64, n_nodes, size(e.coordinates, 2))

#       nf  = zeros(Float64, n_nodes, size(x, 1))
#       nfx = zeros(Float64, n_nodes, size(x, 1))
#       nfy = zeros(Float64, n_nodes, size(x, 1))

#       ReferenceFEs.vander2d!(A, Ax, Ay, e.coordinates, degree)
#       ReferenceFEs.vander2d!(nf, nfx, nfy, x, degree)

#       Ns = (A \ nf)'

#       fn = polyval2d.(e.coordinates[1, :], e.coordinates[2, :], (poly_coeffs,))
      
#       finterpolated = dot(Ns[1, :], fn)
#       @test expected ≈ finterpolated
#     end
#   end

#   @testset "Test SimplexTri grad interpolation" begin
#     x = generate_random_points_in_triangle(1)
#     for degree in 1:max_degree
#       e = ReferenceFEStencil(SimplexTri(), degree)
#       n_nodes = (degree + 1) * (degree + 2) ÷ 2

#       # poly_coeffs = UpperTriangular(ones(Float64, degree + 1, degree + 1))
#       poly_coeffs = UpperTriangular(ones(Float64, degree + 1, degree + 1))
#       poly_coeffs = map(x -> reverse(x), eachrow(poly_coeffs))
#       poly_coeffs = mapreduce(permutedims, vcat, poly_coeffs)
#       expected = polyval2d(x[1], x[2], poly_coeffs)

#       # get shape functions
#       A  = zeros(Float64, n_nodes, size(e.coordinates, 2))
#       Ax = zeros(Float64, n_nodes, size(e.coordinates, 2))
#       Ay = zeros(Float64, n_nodes, size(e.coordinates, 2))

#       nf  = zeros(Float64, n_nodes, size(x, 1))
#       nfx = zeros(Float64, n_nodes, size(x, 1))
#       nfy = zeros(Float64, n_nodes, size(x, 1))

#       ReferenceFEs.vander2d!(A, Ax, Ay, e.coordinates, degree)
#       ReferenceFEs.vander2d!(nf, nfx, nfy, x, degree)

#       dshapes_x = (A \ nfx)'
#       dshapes_y = (A \ nfy)'

#       expected_dx = dpolyval2d(x[1], x[2], poly_coeffs, 1)
#       expected_dy = dpolyval2d(x[1], x[2], poly_coeffs, 2)

#       fn = polyval2d.(e.coordinates[1, :], e.coordinates[2, :], (poly_coeffs,))

#       dfinterpolated_x = dot(dshapes_x[1, :], fn)
#       dfinterpolated_y = dot(dshapes_y[1, :], fn)

#       @test expected_dx ≈ dfinterpolated_x
#       @test expected_dy ≈ dfinterpolated_y
#     end
#   end
# end

# @testset ExtendedTestSet "Quadratures.jl - SimplexTri implementation" begin
#   @testset "Test quad weight positivity" begin
#     for degree in [1, 2, 3, 4, 5, 6]
#       q_rule = Quadrature(SimplexTri(), degree)
#       for w in q_rule.ws
#         @test w > 0.
#       end
#     end
#   end

#   @testset "Test quad point in domain" begin
#     for degree in [1, 2, 3, 4, 5, 6]
#       q_rule = Quadrature(SimplexTri(), degree)
#       for ξ in eachcol(q_rule.ξs)
#         @test is_inside_triangle(ξ)
#       end
#     end
#   end

#   @testset "Test triangle quadrature exactness" begin
#     for degree in [1, 2, 3, 4, 5, 6]
#       q_rule = Quadrature(SimplexTri(), degree)
#       for i in 1:degree
#         for j in 1:degree - i
#           quad_answer = map((ξ, w) -> w * ξ[1]^i * ξ[2]^j, eachcol(q_rule.ξs), q_rule.ws) |> sum
#           exact = integrate_2d_monomial_on_triangle(i, j)
#           @test quad_answer ≈ exact
#         end
#       end
#     end
#   end
# end
