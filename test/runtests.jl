# using Adapt
# using AMDGPU
using Aqua
# using CUDA
using Distributions
using Exodus
using LinearAlgebra
using ReferenceFiniteElements
using StaticArrays
# using Symbolics
using Test

function is_inside_element(::ReferenceFiniteElements.AbstractEdge, point)
  return (point[1] >= -1.) && (point[1] <= 1.)
end

function is_inside_surface_element(::ReferenceFiniteElements.AbstractEdge, point)
  return (point[1] >= 0.) && (point[1] <= 1.)
end

# function is_inside_element(::ReferenceFiniteElements.AbstractHex, point)
#   return (point[1] >= -1.) && (point[1] <= 1.) &&
#          (point[2] >= -1.) && (point[2] <= 1.) &&
#          (point[3] >= -1.) && (point[3] <= 1.)
# end

function is_inside_element(::ReferenceFiniteElements.AbstractQuad, point)
  return (point[1] >= -1.) && (point[1] <= 1.) &&
         (point[2] >= -1.) && (point[2] <= 1.)
end

function is_inside_element(::ReferenceFiniteElements.AbstractTri, point)
  return (point[1] >= 0.) && (point[1] <= 1.) &&
         (point[2] >= 0.) && (point[2] <= 1. - point[1])
end

# function is_inside_element(::ReferenceFiniteElements.AbstractTet, point)
#   return (point[1] >= 0.) && (point[1] <= 1.) &&
#          (point[2] >= 0.) && (point[2] <= 1. - point[1]) && 
#          (point[3] >= 0.) && (point[3] <= 1. - point[2])
# end

q_weight_sum(::ReferenceFiniteElements.AbstractVertex) = 1.
q_weight_sum(e::ReferenceFiniteElements.AbstractEdge) = e.shifted ? 1. : 2. 
# q_weight_sum(::ReferenceFiniteElements.AbstractHex) = 8.
q_weight_sum(::ReferenceFiniteElements.AbstractQuad) = 4.
# q_weight_sum(::ReferenceFiniteElements.AbstractTet) = 1. / 6.
q_weight_sum(::ReferenceFiniteElements.AbstractTri) = 0.5
# q_weight_sum(::ReferenceFiniteElements.AbstractVertex) = 1.

function test_q_points_inside_element(re::ReferenceFE) 
  for n in 1:num_cell_quadrature_points(re)
    @test is_inside_element(re.element, cell_quadrature_point(re, n))
  end

  # TODO re-enable
  # @test all(is_inside_element.((re.element,), surface_quadrature_points(re)))
  # for f in 1:num_faces(re.element)
  #   for n in 1:num_quadrature_points(surface_element(re.element))
  #     if typeof(re.element) <: ReferenceFiniteElements.AbstractTri
  #       @test is_inside_element(ReferenceFiniteElements.surface_element(re.element), 2. * surface_quadrature_point(re, n, f) .- 1.)
  #     else
  #       @test is_inside_element(ReferenceFiniteElements.surface_element(re.element), surface_quadrature_point(re, n, f))
  #     end
  #   end
  # end
end

function test_q_weight_positivity(re::ReferenceFE) 
  for n in 1:num_cell_quadrature_points(re)
    @test cell_quadrature_weight(re, n) > 0.
  end
  for f in 1:num_boundaries(re)
    for n in 1:num_surface_quadrature_points(re)
      @test surface_quadrature_weight(re, n, f) > 0.
    end
  end
end

function test_q_weight_sum(re::ReferenceFE) 
  qwt_sum = 0.
  for q in 1:num_cell_quadrature_points(re)
    qwt_sum += cell_quadrature_weight(re, q)
  end
  @test q_weight_sum(re.element) ≈ qwt_sum
  for f in 1:num_boundaries(re)
    qwt_sum = 0.
    for q in 1:num_surface_quadrature_points(re)
      qwt_sum += surface_quadrature_weight(re, q, f)
    end

    @test q_weight_sum(boundary_element(re)) ≈ qwt_sum
  end
end 

# TODO need a test on exact integration for quadrature rules
function test_partition_of_unity_on_values(re::ReferenceFE) 
  for n in 1:num_cell_quadrature_points(re)
    N = cell_shape_function_value(re, n)
    @test sum(N) ≈ one(eltype(N))
  end

  # TODO re-enable 
  # for f in 1:num_facets(re)
  #   for n in 1:num_surface_quadrature_points(re)
  #     N = surface_shape_function_value(re, n, f)
  #     @test sum(N) ≈ one(eltype(N))
  #   end
  # end
end

function test_partition_of_unity_on_gradients(re::ReferenceFE) 
  for n in 1:num_cell_quadrature_points(re)
    ∇N = cell_shape_function_gradient(re, n)
    @test isapprox(sum(∇N), 0, atol=5e-12)
  end

  for f in 1:num_boundaries(re)
    for n in 1:num_surface_quadrature_points(re)
      ∇N = surface_shape_function_gradient(re, n, f)
      @test isapprox(sum(∇N), 0, atol=5e-11)
    end
  end
end

function test_partition_of_unity_on_hessians(re::ReferenceFE)
  for n in 1:num_cell_quadrature_points(re)
    ∇∇N = cell_shape_function_hessian(re, n)
    @test isapprox(sum(∇∇N), 0, atol=5e-10)
  end

  for f in 1:num_boundaries(re.element)
    for n in 1:num_surface_quadrature_points(re)
      ∇∇N = surface_shape_function_hessian(re, n, f)
      @test isapprox(sum(∇∇N), 0, atol=5e-10)
    end
  end
end

function test_kronecker_delta_property(re::ReferenceFE) 
  Xs = dof_coordinates(re)
  if polynomial_degree(re) == 0
    @test isapprox(hcat(map(x -> 
      ReferenceFiniteElements.shape_function_value(re.element, Xs, x), eachcol(Xs))...), 
      fill(1., (size(Xs, 2), length(Xs))), 
      atol=5e-14
    )
  else
    if dimension(re) == 1
      @test isapprox(hcat(map(x -> 
        ReferenceFiniteElements.shape_function_value(re.element, Xs, x[1]), eachcol(Xs))...), 
        I, 
        atol=5e-12
      )
    else
      @test isapprox(hcat(map(x -> 
        ReferenceFiniteElements.shape_function_value(re.element, Xs, x), eachcol(Xs))...), 
        I, 
        atol=5e-12
      )
    end
  end
  # TODO add test on surface shape functions
  # will need to first calculate shape function then 
  # index based on the edges, since non-edge dofs will be zero 
end

function test_topology_interface_vertex()
  re = Vertex()
  @test boundary_element(re) === nothing
  @test boundary_normals(re) ≈ Matrix{Float64}(undef, 3, 0)
  @test cell_vertices(re) == [1]
  @test dimension(re) == 0
  @test edge_vertices(re) == Matrix{Int}(undef, 0, 0)
  @test face_vertices(re) == Matrix{Int}(undef, 0, 0)
  @test num_boundaries(re) == 0
  @test num_edges(re) == 0
  @test num_faces(re) == 0
  @test num_vertices_per_cell(re) == 1
  @test num_vertices_per_edge(re) == 0
  @test num_vertices_per_face(re) == 0
  @test vertex_coordinates(re) ≈ [0. 0. 0.]' |> collect
end

function test_topology_interface_edge(interp_type, p, shifted)
  re = Edge{interp_type, p}(; shifted = shifted)
  @test boundary_element(re) == Vertex()
  @test boundary_normals(re) == [-1. 0. 0.; 1. 0. 0.]' |> collect
  @test cell_vertices(re) == [1, 2]
  @test dimension(re) == 1
  @test edge_vertices(re) == [1 2]' |> collect
  @test face_vertices(re) == Matrix{Int}(undef, 0, 0)
  @test num_boundaries(re) == 2
  @test num_edges(re) == 1
  @test num_faces(re) == 0
  @test num_vertices_per_cell(re) == 2
  @test num_vertices_per_edge(re) == 2
  @test num_vertices_per_face(re) == 0

  if re.shifted
    @test vertex_coordinates(re) ≈ [0. 0. 0; 1. 0. 0.]'
  else
    @test vertex_coordinates(re) ≈ [-1. 0. 0.; 1. 0. 0.]' |> collect
  end
end

function test_topology_interface_quad(interp_type, p)
  re = Quad{interp_type, p}()
  @test boundary_element(re) == Edge{interp_type, p}()
  @test boundary_normals(re) ≈ [
     0. 1. 0. -1.;
    -1. 0. 1.  0.;
     0. 0. 0.  0.
  ]
  @test cell_vertices(re) == 1:4 |> collect
  @test dimension(re) == 2
  @test edge_vertices(re) == [
    1 2 3 4;
    2 3 4 1
  ]
  @test face_vertices(re) == 1:4 |> collect
  @test num_boundaries(re) == 4
  @test num_edges(re) == 4
  @test num_faces(re) == 1
  @test num_vertices_per_cell(re) == 4
  @test num_vertices_per_edge(re) == 2
  @test num_vertices_per_face(re) == 4
  @test vertex_coordinates(re) ≈ [
    -1.  1. 1. -1.;
    -1. -1. 1.  1.;
     0.  0. 0.  0.
  ]
end

function test_topology_interface_tri(interp_type, p)
  re = Tri{interp_type, p}()
  @test boundary_element(re) == Edge{interp_type, p}(; shifted = true)
  @test boundary_normals(re) ≈ [
     0. 1. / sqrt(2.) -1.;
    -1. 1. / sqrt(2.)  0.;
     0. 0.             0.
  ]
  @test cell_vertices(re) == 1:3 |> collect
  @test dimension(re) == 2
  @test edge_vertices(re) == [
    1 2 3;
    2 3 1
  ]
  @test face_vertices(re) == 1:3 |> collect
  @test num_boundaries(re) == 3
  @test num_edges(re) == 3
  @test num_faces(re) == 1
  @test num_vertices_per_cell(re) == 3
  @test num_vertices_per_edge(re) == 2
  @test num_vertices_per_face(re) == 3
  @test vertex_coordinates(re) ≈ [
    0. 1. 0.;
    0. 0. 1.;
    0. 0. 0.
  ]
  # ev = edge_vertices(fe)
  # @test ev[1:2, 1] == [1, 2]
  # @test ev[1:2, 2] == [2, 3]
  # @test ev[1:2, 3] == [3, 1]
  # if p > 2
  #   offset = 4
  #   for n in 1:3
  #     @test ev[3:end, n] == offset:offset + p - 2
  #     offset += p - 1
  #   end
  # end
  # @test face_vertices(fe) == 1:num_vertices(fe) |> collect
  # if interp_type == Lagrange
  #   if p == 0
  #     @test interior_vertices(fe) == Int[1]
  #   elseif p == 1
  #     @test interior_vertices(fe) == Int[]
  #   else
  #     offset = 3 + 3 * (p - 1) + 1
  #     @test interior_vertices(fe) == offset:offset + (p - 1) * (p - 2) / 2 - 1 |> collect
  #   end
  # end
  # @test num_edges(fe) == 3
  # @test num_faces(fe) == 1
  # if interp_type == Lagrange
  #   @test num_interior_vertices(fe) == (p - 1) * (p - 2) / 2
  # end
  # if interp_type == Lagrange
  #   @test num_vertices(fe) == (p + 1) * (p + 2) / 2
  # end
  # if interp_type == Lagrange
  #   @test num_vertices_per_edge(fe) == p + 1
  # end
  # @test num_vertices_per_face(fe) == (p + 1) * (p + 2) / 2
  # @test surface_element(fe) == Edge{interp_type, p}(; shifted = surface_element(fe).shifted)
  # @test surface_normals(fe) ≈ [
  #    0. 1. / sqrt(2.) -1.;
  #   -1. 1. / sqrt(2.)  0.;
  #    0. 0.             0.
  # ]
  # # TODO vertex coordinates
  # if interp_type == Lagrange
  #   if p == 0
  #     @test vertex_coordinates(fe) ≈ zeros(3, 1)
  #   else
  #     coords = vertex_coordinates(fe)
  #     edge_coords = vertex_coordinates(surface_element(fe))
  #     @test coords[:, 1] ≈ [0., 0., 0.]
  #     @test coords[:, 2] ≈ [1., 0., 0.]
  #     @test coords[:, 3] ≈ [0., 1., 0.]

  #     if p > 1
  #       # test faces
  #       offset = 4
  #       # face 1
  #       for n in 1:p - 1
  #         @test coords[:, offset + n - 1] ≈ [edge_coords[1, n + 2], -1., 0.]
  #       end
  #       offset += p - 1
  #       # face 2
  #       for n in 1:p - 1
  #         @test coords[:, offset + n - 1] ≈ [edge_coords[1, n + 2], 1. - edge_coords[1, n + 2], 0.]
  #       end
  #       offset += p - 1
  #       # face 3
  #       for n in 1:p - 1
  #         @test coords[:, offset + n - 1] ≈ [-1., edge_coords[1, n + 2], 0.]
  #       end
  #       offset += p - 1

  #       # TODO fix
  #       # test interiors
  #       offset = 3 + 3 * (p - 1) + 1
  #       carry = 1
  #       for n in 1:p - 1
  #         for m in 1:p - 1 - n
  #           @test coords[:, offset + carry - 1] ≈ [edge_coords[1, m + 2], edge_coords[1, n + 2], 0.]
  #           carry += 1
  #         end
  #       end
  #     end
  #   end
  # end
end

function test_topology_interface_hex(interp_type, p)
  re = Hex{interp_type, p}()
  @test boundary_element(re) == Quad{interp_type, p}()
  @test boundary_normals(re) ≈ [
    0. 1. 0. -1.  0. 0.
    0. 0. 0.  0. -1. 1.
   -1. 0. 1.  0.  0. 0.
  ]
  @test cell_vertices(re) == 1:8 |> collect
  @test dimension(re) == 3
  @test edge_vertices(re) == [
    1 2 3 4 5 6 7 8 1 2 3 4
    2 3 4 1 6 7 8 5 5 6 7 8
  ]
  @test face_vertices(re) == [
    1 2 3 1 1 5
    2 3 4 5 4 6
    6 7 8 8 3 7
    5 6 7 4 2 8
  ]
  @test num_boundaries(re) == 6
  @test num_edges(re) == 12
  @test num_faces(re) == 6
  @test num_vertices_per_cell(re) == 8
  @test num_vertices_per_edge(re) == 2
  @test num_vertices_per_face(re) == 4
  @test vertex_coordinates(re) ≈ [
    -1.  1.  1. -1. -1.  1. 1. -1.
    -1. -1.  1.  1. -1. -1. 1.  1.
    -1. -1. -1. -1.  1.  1. 1.  1.
  ]
end

function test_topology_interface_tet(interp_type, p)
  re = Tet{interp_type, p}()
  @test boundary_element(re) == Tri{interp_type, p}()
  @test boundary_normals(re) ≈ [
    0. 1. / sqrt(3.) -1.  0.
   -1. 1. / sqrt(3.)  0.  0.
    0. 1. / sqrt(3.)  0. -1.
 ]
  @test cell_vertices(re) == 1:4 |> collect
  @test dimension(re) == 3
  @test edge_vertices(re) == [
    1 2 3 1 2 3
    2 3 1 4 4 4
  ]
  @test face_vertices(re) == [
    1 2 1 1
    2 3 4 3
    4 4 3 2
  ]
  @test num_boundaries(re) == 4
  @test num_edges(re) == 6
  @test num_faces(re) == 4
  @test num_vertices_per_cell(re) == 4
  @test num_vertices_per_edge(re) == 2
  @test num_vertices_per_face(re) == 3
  @test vertex_coordinates(re) ≈ [
    0. 1. 0. 0.
    0. 0. 1. 0.
    0. 0. 0. 1.
  ]
end

function test_dof_interface_vertex()
  re = Vertex()
  @test boundary_dofs(re) == Matrix{Int}(undef, 0, 0)
  @test dof_coordinates(re) ≈ [0. 0. 0.]' |> collect
  @test num_cell_dofs(re) == 1
  @test num_interior_dofs(re) == 0
end

function test_dof_interface_edge(interp_type::Type{Lagrange}, p, shifted)
  re = Edge{interp_type, p}(; shifted = shifted)
  if re.shifted
    x_min = 0.
  else
    x_min = -1.
  end

  @test boundary_dofs(re) == [1 2]' |> collect
  coords = dof_coordinates(re)
  if p == 0
    @test coords ≈ zeros(1, 1)
  elseif p == 1
    @test coords ≈ [x_min 1.]
  else
    @test coords[1] ≈ x_min
    @test coords[2] ≈ 1.
    coords_rng = LinRange(x_min, 1., p + 1)
    for n in 1:p - 2
      @test coords[n + 2] ≈ coords_rng[n + 1]
    end
  end
  if p < 2
    @test interior_dofs(re) == Int[]
  else
    @test interior_dofs(re) == 3:p + 1
  end
  @test num_cell_dofs(re) == p + 1
  if p < 2
    @test num_interior_dofs(re) == 0
  else
    @test num_interior_dofs(re) == p - 1
  end
end

function test_dof_interface_quad(interp_type::Type{Lagrange}, p)
  re = Quad{interp_type, p}()
  bd = boundary_dofs(re)
  @test bd[1:2, 1] == [1, 2]
  @test bd[1:2, 2] == [2, 3]
  @test bd[1:2, 3] == [3, 4]
  @test bd[1:2, 4] == [4, 1]
  if p > 2
    offset = 5
    for n in 1:4
      @test bd[3:end, n] == offset:offset + p - 2
      offset += p - 1
    end
  end
  coords = dof_coordinates(re)
  if p == 0
    @test coords ≈ zeros(2, 1)
  elseif p == 1
    @test coords ≈ [
      -1.  1. 1. -1.
      -1. -1. 1.  1.
    ]
  else
    edge_coords = dof_coordinates(boundary_element(re))
    # test faces
    offset = 5
    # face 1
    for n in 1:p - 1
      @test coords[:, offset + n - 1] ≈ [edge_coords[1, n + 2], -1.]
    end
    offset += p - 1
    # face 2
    for n in 1:p - 1
      @test coords[:, offset + n - 1] ≈ [1., edge_coords[1, n + 2]]
    end
    offset += p - 1
    # face 3
    for n in 1:p - 1
      @test coords[:, offset + n - 1] ≈ [edge_coords[1, n + 2], 1.]
    end
    offset += p - 1
    # face 4
    for n in 1:p - 1
      @test coords[:, offset + n - 1] ≈ [-1., edge_coords[1, n + 2]]
    end
    offset += p - 1

    # test interiors
    offset = 4 + 4 * (p - 1) + 1
    carry = 1
    for n in 1:p - 1
      for m in 1:p - 1
        @test coords[:, offset + carry - 1] ≈ [edge_coords[1, m + 2], edge_coords[1, n + 2]]
        carry += 1
      end
    end
  end
  if p < 2
    @test interior_dofs(re) == Int[]
  else
    @test interior_dofs(re) == 4 + 4 * (p - 1) + 1:4 + 4 * (p - 1) + (p - 1) * (p - 1) |> collect
  end
  @test num_cell_dofs(re) == (p + 1) * (p + 1)
  @test num_interior_dofs(re) == (p - 1) * (p - 1)
end

# TODO finish this up by testing boundary_dofs and interior_dofs
function test_dof_interface_tri(interp_type::Type{Lagrange}, p)
  re = Tri{interp_type, p}()
  coords = dof_coordinates(re)
  if p == 0
    @test coords ≈ zeros(2, 1)
  elseif p == 1
    @test coords ≈ [
      0. 1. 0.;
      0. 0. 1.
    ]
  else
    edge_coords = dof_coordinates(boundary_element(re))
    # test faces
    offset = 4
    # face 1
    for n in 1:p - 1
      @test coords[:, offset + n - 1] ≈ [edge_coords[1, n + 2], -1.]
    end
    offset += p - 1
    # face 2
    for n in 1:p - 1
      @test coords[:, offset + n - 1] ≈ [edge_coords[1, n + 2], 1. - edge_coords[1, n + 2]]
    end
    offset += p - 1
    # face 3
    for n in 1:p - 1
      @test coords[:, offset + n - 1] ≈ [-1., edge_coords[1, n + 2]]
    end
    offset += p - 1

    # TODO fix this up
    # test interiors
    offset = 3 + 3 * (p - 1) + 1
    carry = 1
    for n in 1:p - 1
      for m in 1:p - 1 - n
        @test coords[:, offset + carry - 1] ≈ [edge_coords[1, m + 2], edge_coords[1, n + 2]]
        carry += 1
      end
    end
  end
  @test num_cell_dofs(re) == (p + 1) * (p + 2) ÷ 2
  if p < 3
    @test num_interior_dofs(re) == 0
  else
    @test num_interior_dofs(re) == (p - 1) * (p - 2) ÷ 2
  end
end

function test_ref_fe(el_type, interp_type, p, q_rule; kwargs...)
  el = el_type{interp_type, p}(; kwargs...)
  @testset "$el - $q_rule Tests" begin
    re = ReferenceFE(el, q_rule)
    test_q_points_inside_element(re)
    # TODO eventually disable this for tets
    # if el_type <: ReferenceFiniteElements.AbstractTet || q_order > 1
    #   # do nothing here
    # else
      test_q_weight_positivity(re)
    # end
    test_q_weight_sum(re)
    test_partition_of_unity_on_values(re)
    test_partition_of_unity_on_gradients(re)
    test_partition_of_unity_on_hessians(re)
    # test_kronecker_delta_property(re)
  end
end

# Topology interface tests
@testset "Topology interface - Vertex" begin
  test_topology_interface_vertex()
end

@testset "Topology interface - Edge" begin
  test_topology_interface_edge(Lagrange, 0, false)
  test_topology_interface_edge(Lagrange, 0, true)
end

@testset "Topology interface - Quad" begin
  test_topology_interface_quad(Lagrange, 0)
end

@testset "Topology interface - Tri" begin
  test_topology_interface_tri(Lagrange, 0)
end

@testset "Topology interface - Hex" begin
  test_topology_interface_hex(Lagrange, 0)
end

@testset "Topology interface - Tet" begin
  test_topology_interface_tet(Lagrange, 0)
end

# Dof interface tests
@testset "Dof Interface - Vertex" begin
  test_dof_interface_vertex()
end

@testset "Dof Interface - Edge" begin
  for p in 0:5
    test_dof_interface_edge(Lagrange, p, false)
    test_dof_interface_edge(Lagrange, p, true)
  end
end

@testset "Dof Interface - Quad" begin
  for p in 0:5
    test_dof_interface_quad(Lagrange, p)
  end
end

@testset "Dof Interface - Tri" begin
  for p in 0:5
    test_dof_interface_tri(Lagrange, p)
  end
end

# ReferenceFE tests

# TODO add Vertex tests

for p in 1:5
  q_rule = GaussLobattoLegendre(p)
  test_ref_fe(Edge, Lagrange, p, q_rule; shifted = false)
  test_ref_fe(Edge, Lagrange, p, q_rule; shifted = true)
end

for p in 1:5
  q_rule = GaussLobattoLegendre(p)
  test_ref_fe(Quad, Lagrange, p, q_rule)
end

for p in 1:5
  q_rule = GaussLobattoLegendre(p)
  test_ref_fe(Tri, Lagrange, p, q_rule)
end

# end

# @testset "ReferenceFEs" begin
#   test_ref_fes()
# end

# extension tests
# include("TestAdaptExt.jl")
include("TestExodusExt.jl")
# include("TestSymbolicsExt.jl")

@testset "ExodusExt" begin
  test_exodus_ext()
end

# @testset "SymbolicExt" begin
#   test_symbolic_fe_fes()
# end

# @testset "Aqua Tests" begin
#   Aqua.test_all(ReferenceFiniteElements; ambiguities=false)
# end
