# abstract methods
surface_element(::AbstractTri{V, I, P, Q}) where {V, I, P, Q} = Edge{I, P, Q}()

function element_edge_vertices(e::AbstractTri, backend)
  if polynomial_degree(e) == 0
    n_verts_per_edge = 2
  else
    n_verts_per_edge = num_vertices_per_edge(e)
  end

  # edges = Matrix{Int64}(undef, num_vertices_per_edge(e), 3)
  edges = Matrix{Int64}(undef, n_verts_per_edge, 3)
  edge_nodes = 1:polynomial_degree(e) + 1

  # corners
  edges[1, 1] = 1
  edges[2, 1] = 2
  edges[1, 2] = 2
  edges[2, 2] = 3
  edges[1, 3] = 3
  edges[2, 3] = 1

  # middle edge guys
  curr_node = 4
  if length(edge_nodes) > 2
    # for n in 2:length(edge_nodes) - 1
    for n in 3:length(edge_nodes)
      for m in 1:3
        edges[n, m] = curr_node
        curr_node = curr_node + 1
      end
    end
  end

  return map(x -> convert_to(backend, x...), eachcol(edges))
end

function element_face_vertices(e::AbstractTri, backend::ArrayBackend)
  temp = convert_to(backend, zeros(Int, 0)...)
  face_nodes = [temp]
  return face_nodes
end

# TODO
function element_interior_nodes(e::AbstractTri, backend::ArrayBackend)
  nodes = Vector{Int64}(undef, num_interior_vertices(e))
  # for n in axes(nodes, 1)
  #   nodes[n] = 4 * num_nodes_per_edge(e) - 4 + n
  # end
  return nodes
end

function nodal_coordinates(e::AbstractTri, backend)
  edge_coords = nodal_coordinates(surface_element(e), backend)
  edge_coords = map(x -> 0.5 * (x .+ 1), edge_coords)
  coords = Matrix{Float64}(undef, 2, num_vertices(e))

  # corner nodes
  coords[1, 1] = 0.
  coords[2, 1] = 0.
  coords[1, 2] = 1.
  coords[2, 2] = 0.
  coords[1, 3] = 0.
  coords[2, 3] = 1. 

  # edge mid point nodes
  curr_node = 4
  @views for val in edge_coords[2:end - 1]
    coords[1, curr_node] = val[1]
    coords[2, curr_node] = 0.
    curr_node = curr_node + 1

    coords[1, curr_node] = val[1]
    coords[2, curr_node] = 1. - val[1]
    curr_node = curr_node + 1

    coords[1, curr_node] = 0.
    coords[2, curr_node] = val[1]
    curr_node = curr_node + 1
  end

  # face interior nodes
  for val in Iterators.product(edge_coords[2:end - 1], edge_coords[2:end - 1])
    # TODO need to check this...
    if val[2][1] >= 1. - val[1][1]
      continue
    end 
    coords[1, curr_node] = val[1][1]
    coords[2, curr_node] = val[2][1]
    curr_node = curr_node + 1
  end

  # return coords
  return map(x -> convert_to_vector_coords(e, backend, x...), eachcol(coords))
end

# TODO this needs to be fixed maybe?
function surface_nodal_coordinates(e::AbstractTri, backend)
  coords = nodal_coordinates(e, backend)
  coords = map(x -> 0.5 * (x .+ 1), coords)
  edges = element_edge_vertices(e, backend)
  surf_coords = mapreduce(x -> coords[x] |> collect, hcat, edges)
  return surf_coords
end

function quadrature_points_and_weights(e::AbstractTri, backend::ArrayBackend)
  if quadrature_degree(e) == 1
    ξs = Matrix{Float64}(undef, 2, 1)
    ξs[:, 1] = [1. / 3., 1. / 3.]
    ws = [0.5]
  elseif quadrature_degree(e) == 2
    ξs = Matrix{Float64}(undef, 2, 3)
    ξs[:, 1] = [2. / 3., 1. / 6.]
    ξs[:, 2] = [1. / 6., 2. / 3.]
    ξs[:, 3] = [1. / 6., 1. / 6.]
    ws    = [1. / 6., 1. / 6., 1. / 6.]
  elseif quadrature_degree(e) <= 4
    ξs = Matrix{Float64}(undef, 2, 6)
    ξs[:, 1] = [1.081030181680700E-01, 4.459484909159650E-01]
    ξs[:, 2] = [4.459484909159650E-01, 1.081030181680700E-01]
    ξs[:, 3] = [4.459484909159650E-01, 4.459484909159650E-01]
    ξs[:, 4] = [8.168475729804590E-01, 9.157621350977100E-02]
    ξs[:, 5] = [9.157621350977100E-02, 8.168475729804590E-01]
    ξs[:, 6] = [9.157621350977100E-02, 9.157621350977100E-02]

    ws = [
      1.116907948390055E-01,
      1.116907948390055E-01,
      1.116907948390055E-01,
      5.497587182766100E-02,
      5.497587182766100E-02,
      5.497587182766100E-02
    ]
  elseif quadrature_degree(e) <= 5
    ξs = Matrix{Float64}(undef, 2, 7)
    ξs[:, 1] = [3.33333333333333E-01, 3.33333333333333E-01]
    ξs[:, 2] = [5.97158717897700E-02, 4.70142064105115E-01]
    ξs[:, 3] = [4.70142064105115E-01, 5.97158717897700E-02]
    ξs[:, 4] = [4.70142064105115E-01, 4.70142064105115E-01]
    ξs[:, 5] = [7.97426985353087E-01, 1.01286507323456E-01]
    ξs[:, 6] = [1.01286507323456E-01, 7.97426985353087E-01]
    ξs[:, 7] = [1.01286507323456E-01, 1.01286507323456E-01]

    ws = [
      1.12500000000000E-01,
      6.61970763942530E-02,
      6.61970763942530E-02,
      6.61970763942530E-02,
      6.29695902724135E-02,
      6.29695902724135E-02,
      6.29695902724135E-02
    ]
  elseif quadrature_degree(e) <= 6
    ξs = Matrix{Float64}(undef, 2, 12)
    ξs[:, 1]  = [5.01426509658179E-01, 2.49286745170910E-01]
    ξs[:, 2]  = [2.49286745170910E-01, 5.01426509658179E-01]
    ξs[:, 3]  = [2.49286745170910E-01, 2.49286745170910E-01]
    ξs[:, 4]  = [8.73821971016996E-01, 6.30890144915020E-02]
    ξs[:, 5]  = [6.30890144915020E-02, 8.73821971016996E-01]
    ξs[:, 6]  = [6.30890144915020E-02, 6.30890144915020E-02]
    ξs[:, 7]  = [5.31450498448170E-02, 3.10352451033784E-01]
    ξs[:, 8]  = [6.36502499121399E-01, 5.31450498448170E-02]
    ξs[:, 9]  = [3.10352451033784E-01, 6.36502499121399E-01]
    ξs[:, 10] = [5.31450498448170E-02, 6.36502499121399E-01]
    ξs[:, 11] = [6.36502499121399E-01, 3.10352451033784E-01]
    ξs[:, 12] = [3.10352451033784E-01, 5.31450498448170E-02]

    ws = [
      5.83931378631895E-02,
      5.83931378631895E-02,
      5.83931378631895E-02,
      2.54224531851035E-02,
      2.54224531851035E-02,
      2.54224531851035E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02
    ]
  end
  return map(x -> convert_to(backend, x...), eachcol(ξs)), ws
end

function surface_quadrature_points_and_weights(e::AbstractTri, backend::ArrayBackend)
  ξs, ws = quadrature_points_and_weights(surface_element(e), backend)
  # ξs = (ξs .+ 1.) / 2.
  ξs = map(x -> (x .+ 1.) / 2, ξs)
  ws = ws / 2

  face_1_ξs = map(x -> vcat(x, 0.), ξs)
  face_2_ξs = map(x -> vcat(x, 1. .- x), ξs)
  face_3_ξs = map(x -> vcat(0., x), ξs)

  return [face_1_ξs, face_2_ξs, face_3_ξs], [ws, ws, ws]
end 

# Specific implementations

# Constant tri
"""
$(TYPEDEF)
"""
struct Tri0{I, Q} <: AbstractTri{3, I, 0, Q}
end

# Lagrange implementation
function shape_function_value(e::Tri0{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_vector(e, backend,
    1.
  )
  return Ns
end

function shape_function_gradient(e::Tri0{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_matrix(e, backend,
    0.,
    #
    0.
  )
  return Ns
end

function shape_function_hessian(e::Tri0{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_3d_array(e, backend,
    0., 0.,
    #
    0., 0.
  )
  return Ns
end

# Linear tri
"""
$(TYPEDEF)
"""
struct Tri3{I, Q} <: AbstractTri{3, I, 1, Q}
end

# Lagrange implementation
function shape_function_value(e::Tri3{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_vector(e, backend,
    1. - ξ[1] - ξ[2],
    ξ[1],
    ξ[2]
  )
  return Ns
end

function shape_function_gradient(e::Tri3{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_matrix(e, backend,
    -1., 
    1., 
    0.,
    #
    -1., 
    0., 
    1.
  )
  return Ns
end

function shape_function_hessian(e::Tri3{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_3d_array(e, backend,
    0., 0.,
    0., 0.,
    0., 0.,
    #
    0., 0.,
    0., 0.,
    0., 0.
  )
  return Ns
end

# Quadratic tri
"""
$(TYPEDEF)
"""
struct Tri6{I, Q} <: AbstractTri{6, I, 2, Q}
end

# Lagrange implementation
function shape_function_value(e::Tri6{Lagrange}, X, ξ, backend::ArrayBackend)
  λ = 1. - ξ[1] - ξ[2]
  Ns = convert_to_vector(e, backend,
    λ * (2. * λ - 1.),
    ξ[1] * (2. * ξ[1] - 1.),
    ξ[2] * (2. * ξ[2] - 1.),
    4. * ξ[1] * λ,
    4. * ξ[1] * ξ[2],
    4. * ξ[2] * λ
  )
  return Ns
end

function shape_function_gradient(e::Tri6{Lagrange}, X, ξ, backend::ArrayBackend)
  λ = 1. - ξ[1] - ξ[2]
  Ns = convert_to_matrix(e, backend,
    -1. * (2. * λ - 1.) - 2. * λ,
    (2. * ξ[1] - 1.) + 2. * ξ[1],
    0.0,
    4. * λ - 4. * ξ[1],
    4. * ξ[2],
  -4. * ξ[2],
  #
  -1. * (2. * λ - 1.) - 2. * λ,
    0.0,
    (2. * ξ[2] - 1.) + 2. * ξ[2],
  -4. * ξ[1],
    4. * ξ[1],
    4. * λ - 4. * ξ[2]
  )
  return Ns
end

function shape_function_hessian(e::Tri6{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_3d_array(e, backend,
    4., 4., 
    0., -8., 
    0., 0.,
    4., 0., 
    0., -4., 
    4., -4.,
    #
    4., 0., 
    0., -4., 
    4., -4.,
    4., 0., 
    4.,  0., 
    0., -8.
  )
  return Ns
end

# General Tri
"""
$(TYPEDEF)
"""
struct Tri{V, I, P, Q} <: AbstractTri{V, I, P, Q}
end

# Lagrange implementation
function shape_function_value(e::Tri{V, Lagrange}, X, ξ, backend::ArrayBackend) where V
  coords_x = nodal_coordinates(surface_element(e), backend)
  coords_y = nodal_coordinates(surface_element(e), backend)
  coords_x = map(x -> 0.5 * (x .+ 1), coords_x)
  coords_y = map(x -> 0.5 * (x .+ 1), coords_y)
  coords_t = map((x, y) -> 1. .- x .- y, coords_x, coords_y)
  N_x = shape_function_value(surface_element(e), coords_x, ξ[1], backend)
  N_y = shape_function_value(surface_element(e), coords_y, ξ[2], backend)
  N_t = shape_function_value(surface_element(e), coords_t, 1. - ξ[1] - ξ[2], backend)

  N = Vector{eltype(coords_x[1])}(undef, num_vertices(e))

  # Ns = convert_to_vector(e, backend,
    
  # )

  # corner nodes first
  # N[1] = N_x[1] * N_y[1] * N_t[1]
  # N[2] = N_x[end] * N_y[1] * N_t[1]
  # N[3] = N_x[1] * N_y[1] * N_t[end]

  # # edge nodes next
  # for n in 2:num_nodes_per_edge(e) - 1
  #   k = 3 * (n - 2)
  #   N[3 + k + 1] = N_x[n] * N_y[1] * N_t[1]
  #   N[3 + k + 2] = N_x[end] * N_y[1] * N_t[n]
  #   N[3 + k + 3] = N_x[end] * N_y[n] * N_t[end]
  # end

  # for n in axes(N, 1)
  #   N_r = 2 * ξ[1] / 
  # end

  # interior nodes next

  return convert_to_vector(e, backend, N...)
end

# function shape_function_gradient(e::Tri{Lagrange}, X, ξ, backend::ArrayBackend)
#   Ns = convert_to_matrix(e, backend,
#     -1., 
#     1., 
#     0.,
#     #
#     -1., 
#     0., 
#     1.
#   )
#   return Ns
# end

# function shape_function_hessian(e::Tri{Lagrange}, X, ξ, backend::ArrayBackend)
#   Ns = convert_to_3d_array(e, backend,
#     0., 0.,
#     0., 0.,
#     0., 0.,
#     #
#     0., 0.,
#     0., 0.,
#     0., 0.
#   )
#   return Ns
# end
