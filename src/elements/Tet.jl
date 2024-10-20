# abstract methods
# surface_element(::AbstractTet{V, I, P, Q}) where {V, I, P, Q} = Tri{3, I, P, Q}() # fix this! TODO

# TODO
function element_edge_nodes(e::AbstractTet, backend)
  edges = Matrix{Int64}(undef, num_vertices_per_edge(e), 6)
  edge_nodes = 1:polynomial_degree(e) + 1
  # TODO finish this
  # corners
  # edges[1, 1] = 1
  # edges[end, 1] = 2
  # edges[1, 2] = 2
  # edges[end, 2] = 3
  # edges[1, 3] = 3
  # edges[end, 3] = 1

  # # middle edge guys
  # curr_node = 4
  # if length(edge_nodes) > 2
  #   for n in 2:length(edge_nodes) - 1
  #     for m in 1:3
  #       edges[n, m] = curr_node
  #       curr_node = curr_node + 1
  #     end
  #   end
  # end

  return map(x -> convert_to(backend, x...), eachcol(edges))
end

function element_face_nodes(e::AbstractTet, backend::ArrayBackend)
  temp = convert_to(backend, zeros(Int, 0)...)
  face_nodes = [temp]
  return face_nodes
end

# TODO
function element_interior_nodes(e::AbstractTet, backend::ArrayBackend)
  nodes = Vector{Int64}(undef, num_interior_vertices(e))
  # for n in axes(nodes, 1)
  #   nodes[n] = 4 * num_nodes_per_edge(e) - 4 + n
  # end
  return nodes
end

function nodal_coordinates(e::AbstractTet, backend::ArrayBackend)
  edge_coords = nodal_coordinates(surface_element(surface_element(e)), backend)
  edge_coords = map(x -> 0.5 * (x .+ 1), edge_coords)
  coords = Matrix{Float64}(undef, 3, num_vertices(e))

  # corner nodes
  coords[1, 1] = 0.
  coords[2, 1] = 0.
  coords[3, 1] = 0.
  #
  coords[1, 2] = 1.
  coords[2, 2] = 0.
  coords[3, 2] = 0.
  #
  coords[1, 3] = 0.
  coords[2, 3] = 1.
  coords[3, 3] = 0.
  #
  coords[1, 4] = 0.
  coords[2, 4] = 0.
  coords[3, 4] = 1.
  
  # TODO need to check below
  # edge mid point nodes
  curr_node = 5
  @views for val in edge_coords[2:end - 1]
    coords[1, curr_node] = val[1]
    coords[2, curr_node] = 0.
    coords[3, curr_node] = 0.
    curr_node = curr_node + 1

    coords[1, curr_node] = val[1]
    coords[2, curr_node] = 1. - val[1]
    coords[3, curr_node] = 0.0
    curr_node = curr_node + 1

    coords[1, curr_node] = 0.
    coords[2, curr_node] = val[1]
    coords[3, curr_node] = 0.
    curr_node = curr_node + 1

    coords[1, curr_node] = 0.0
    coords[2, curr_node] = 0.0
    coords[3, curr_node] = val[1]
    curr_node = curr_node + 1

    coords[1, curr_node] = val[1]
    coords[2, curr_node] = 0.0
    coords[3, curr_node] = 1. - val[1]
    curr_node = curr_node + 1

    coords[1, curr_node] = 0.0
    coords[2, curr_node] = val[1]
    coords[3, curr_node] = 1. - val[1]
    curr_node = curr_node + 1
  end

  # cell interior points
  # TODO

  return map(x -> convert_to_vector_coords(e, backend, x...), eachcol(coords))
end

function surface_nodal_coordinates(e::AbstractTet, backend::ArrayBackend)
  # coords = surface_nodal_coordinates(, backend)
  coords = nodal_coordinates(surface_element(e), backend)
  # edges = element_edge_nodes(e, backend)
  face_1_Xs = map(x-> vcat(x[1], 0., x[2]), coords)
  face_2_Xs = map(x-> vcat(x[1], 0., x[2]), coords)
  face_3_Xs = map(x-> vcat(x[1], 0., x[2]), coords)
  face_4_Xs = map(x-> vcat(x[1], 0., x[2]), coords)
  return [face_1_Xs, face_2_Xs, face_3_Xs, face_4_Xs]
end

function quadrature_points_and_weights(e::AbstractTet, backend::ArrayBackend)
  if quadrature_degree(e) == 1
    ξs = Matrix{Float64}(undef, 3, 1)
    ξs[:, 1] = [1. / 4., 1. / 4., 1. / 4.]
    ws = [1. / 6.]
  elseif quadrature_degree(e) == 2
    ξs = Matrix{Float64}(undef, 3, 5)
    ξs[:, 1] = [1. / 4., 1. / 4., 1. / 4.]
    ξs[:, 2] = [1. / 6., 1. / 6., 1. / 6.]
    ξs[:, 3] = [1. / 6., 1. / 6., 1. / 2.]
    ξs[:, 4] = [1. / 6., 1. / 2., 1. / 6.]
    ξs[:, 5] = [1. / 2., 1. / 6., 1. / 6.]

    #
    ws = [
      -2. / 15.
       3. / 40.
       3. / 40.
       3. / 40.
       3. / 40.
    ]
  end
  return map(x -> convert_to(backend, x...), eachcol(ξs)), ws
end

function surface_quadrature_points_and_weights(e::AbstractTet, backend::ArrayBackend)
  ξs, ws = quadrature_points_and_weights(surface_element(e), backend)

  face_1_ξs = map(x -> vcat(x[1], 0., x[2]), ξs)
  face_2_ξs = map(x -> vcat(x[1], x[2], 1. - x[1] - x[2]), ξs)
  face_3_ξs = map(x -> vcat(0., x[1], x[2]), ξs)
  face_4_ξs = map(x -> vcat(x[1], x[2], 0.), ξs)

  return [face_1_ξs, face_2_ξs, face_3_ξs, face_4_ξs], [ws, ws, ws, ws]
end 

# Specific implementations

# Constant tri
"""
$(TYPEDEF)
"""
struct Tet0{I, Q} <: AbstractTet{4, I, 0, Q}
end
surface_element(::Tet0{I, Q}) where {I, Q} = Tri0{I, Q}()

# Lagrange implementation
function shape_function_value(e::Tet0{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_vector(e, backend,
    1.
  )
  return Ns
end

function shape_function_gradient(e::Tet0{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_matrix(e, backend,
    0.,
    #
    0.,
    #
    0.
  )
  return Ns
end

function shape_function_hessian(e::Tet0{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_3d_array(e, backend,
    0., 0., 0.,
    #
    0., 0., 0.,
    #
    0., 0., 0.
  )
  return Ns
end

# Linear Tet
"""
$(TYPEDEF)
"""
struct Tet4{I, Q} <: AbstractTet{4, I, 1, Q}
end
surface_element(::Tet4{I, Q}) where {I, Q} = Tri3{I, Q}()

# Lagrange implementation
function shape_function_value(e::Tet4{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_vector(e, backend,
    1. - ξ[1] - ξ[2] - ξ[3],
    ξ[1],
    ξ[2],
    ξ[3]
  )
  return Ns
end

function shape_function_gradient(e::Tet4{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_matrix(e, backend,
    -1., 
    1., 
    0., 
    0.,
    #
    -1., 
    0., 
    1., 
    0.,
    #
    -1., 
    0., 
    0., 
    1. 
  )
  return Ns
end

function shape_function_hessian(e::Tet4{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_3d_array(e, backend,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    #
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    #
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.
  )
  return Ns
end

# Quadratic Tet
"""
$(TYPEDEF)
"""
struct Tet10{I, Q} <: AbstractTet{10, I, 2, Q}
end
surface_element(::Tet10{I, Q}) where {I, Q} = Tri6{I, Q}()

# Lagrange implementation
function shape_function_value(e::Tet10{Lagrange}, X, ξ, backend::ArrayBackend)
  t0 = 1 - ξ[1] - ξ[2] - ξ[3]
  t1 = ξ[1]
  t2 = ξ[2]
  t3 = ξ[3]
  Ns = convert_to_vector(e, backend,
    t0 * (2 * t0 - 1),
    t1 * (2 * t1 - 1),
    t2 * (2 * t2 - 1),
    t3 * (2 * t3 - 1),
    4 * t0 * t1,
    4 * t1 * t2,
    4 * t2 * t0,
    4 * t0 * t3,
    4 * t1 * t3,
    4 * t2 * t3
  )
  return Ns
end

function shape_function_gradient(e::Tet10{Lagrange}, X, ξ, backend::ArrayBackend)
  t0 = 1 - ξ[1] - ξ[2] - ξ[3]
  t1 = ξ[1]
  t2 = ξ[2]
  t3 = ξ[3]
  Ns = convert_to_matrix(e, backend,
    1-4*t0,
    4*t1-1,
    0,
    0,
    4*(t0-t1),
    4*t2,
  -4*t2,
  -4*t3,
    4*t3,
    0,
  #
    1-4*t0,
    0,
    4*t2-1,
    0,
  -4*t1,
    4*t1,
    4*(t0-t2),
  -4*t3,
    0,
    4*t3,
  #
    1-4*t0,
    0,
    0,
    4*t3-1,
  -4*t1,
    0,
  -4*t2,
    4*(t0-t3),
    4*t1,
    4*t2
  )
  return Ns
end

function shape_function_hessian(e::Tet10{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_3d_array(e, backend,
    4,  4,  0,  0, -8,  0,  0,  0,  0,  0,
    4,  0,  0,  0, -4,  4, -4,  0,  0,  0,
    4,  0,  0,  0, -4,  0,  0, -4,  4,  0,
    4,  0,  0,  0, -4,  4, -4,  0,  0,  0,
    4,  0,  4,  0,  0,  0, -8,  0,  0,  0,
    4,  0,  0,  0,  0,  0, -4, -4,  0,  4,
    4,  0,  0,  0, -4,  0,  0, -4,  4,  0,
    4,  0,  0,  0,  0,  0, -4, -4,  0,  4,
    4,  0,  0,  4,  0,  0,  0, -8,  0,  0,
  )
  return Ns
end
