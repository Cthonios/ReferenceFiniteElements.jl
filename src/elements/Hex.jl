# abstract methods
# surface_element(::AbstractHex{V, I, P, Q}) where {V, I, P, Q} = Quad{Int(cbrt(V)), I, P, Q}()

function element_edge_nodes(e::AbstractHex, backend::ArrayBackend)

end

function element_face_nodes(e::AbstractHex, backend::ArrayBackend)
  faces = Matrix{Int64}(undef, num_vertices_per_face(e), 6)
  edge_nodes = 1:polynomial_degree(e) + 1

  # corners
  # faces[1, 1] = 1
  return faces
end

function element_interior_nodes(e::AbstractHex, backend::ArrayBackend)
  nodes = Vector{Int64}(undef, num_interior_vertices(e))
  for n in axes(nodes, 1)
    nodes[n] = 12 * num_vertices_per_edge(e) - 8 + n
  end
  return nodes
end

function nodal_coordinates(e::AbstractHex, backend::ArrayBackend)
  edge_coords = nodal_coordinates(surface_element(surface_element(e)), backend)
  coords = Matrix{Float64}(undef, 3, num_vertices(e))

  # corner nodes
  coords[1, 1] = -1.
  coords[2, 1] = -1.
  coords[3, 1] = -1.
  #
  coords[1, 2] = 1.
  coords[2, 2] = -1.
  coords[3, 2] = -1.
  #
  coords[1, 3] = 1.
  coords[2, 3] = 1.
  coords[3, 3] = -1.
  #
  coords[1, 4] = -1.
  coords[2, 4] = 1.
  coords[3, 4] = -1.
  #
  coords[1, 5] = -1.
  coords[2, 5] = -1.
  coords[3, 5] = 1.
  #
  coords[1, 6] = 1.
  coords[2, 6] = -1.
  coords[3, 6] = 1.
  #
  coords[1, 7] = 1.
  coords[2, 7] = 1.
  coords[3, 7] = 1.
  #
  coords[1, 8] = -1.
  coords[2, 8] = 1.
  coords[3, 8] = 1.

  # edge mid point nodes
  curr_node = 9

  # face 1
  for (val, val_rev) in zip(edge_coords[2:end - 1], reverse(edge_coords[2:end - 1]))
    coords[1, curr_node] = val[1]
    coords[2, curr_node] = -1.
    coords[3, curr_node] = 1.
    curr_node = curr_node + 1

    coords[1, curr_node] = 1.
    coords[2, curr_node] = val[1]
    coords[3, curr_node] = 1.
    curr_node = curr_node + 1

    coords[1, curr_node] = val_rev[1]
    coords[2, curr_node] = 1.
    coords[3, curr_node] = 1.
    curr_node = curr_node + 1

    coords[1, curr_node] = -1.
    coords[2, curr_node] = val_rev[1]
    coords[3, curr_node] = 1.
    curr_node = curr_node + 1
  end

  # face 2
  for (val, val_rev) in zip(edge_coords[2:end - 1], reverse(edge_coords[2:end - 1]))
    coords[1, curr_node] = 1.
    coords[2, curr_node] = -1.
    coords[3, curr_node] = val[1]
    curr_node = curr_node + 1

    coords[1, curr_node] = 1.
    coords[2, curr_node] = val[1]
    coords[3, curr_node] = -1.
    curr_node = curr_node + 1

    coords[1, curr_node] = 1.
    coords[2, curr_node] = 1.
    coords[3, curr_node] = val_rev[1]
    curr_node = curr_node + 1

    # coords[1, curr_node] = -1.
    # coords[2, curr_node] = val_rev[1]
    # coords[3, curr_node] = -1.
    # curr_node = curr_node + 1
  end

  # face 3
  for (val, val_rev) in zip(edge_coords[2:end - 1], reverse(edge_coords[2:end - 1]))
    coords[1, curr_node] = val_rev[1]
    coords[2, curr_node] = -1.
    coords[3, curr_node] = -1.
    curr_node = curr_node + 1

    coords[1, curr_node] = -1.
    coords[2, curr_node] = val[1]
    coords[3, curr_node] = -1.
    curr_node = curr_node + 1

    coords[1, curr_node] = val[1]
    coords[2, curr_node] = 1.
    coords[3, curr_node] = -1.
    curr_node = curr_node + 1

    # coords[1, curr_node] = 1.
    # coords[2, curr_node] = val_rev[1]
    # coords[3, curr_node] = -1.
    # curr_node = curr_node + 1
  end

  # face 4
  for (val, val_rev) in zip(edge_coords[2:end - 1], reverse(edge_coords[2:end - 1]))
    # coords[1, curr_node] = -1.
    # coords[2, curr_node] = val[1]
    # coords[3, curr_node] = 1.
    # curr_node = curr_node + 1

    coords[1, curr_node] = -1.
    coords[2, curr_node] = 1.
    coords[3, curr_node] = val[1]
    curr_node = curr_node + 1

    # coords[1, curr_node] = -1.
    # coords[2, curr_node] = val_rev[1]
    # coords[3, curr_node] = -1.
    # curr_node = curr_node + 1

    coords[1, curr_node] = -1.
    coords[2, curr_node] = -1.
    coords[3, curr_node] = val_rev[1]
    curr_node = curr_node + 1
  end

  # now for interior face points
  # face 1
  for val in Iterators.product(edge_coords[2:end - 1], edge_coords[2:end - 1])
    # face 1
    coords[1, curr_node] = val[1][1]
    coords[2, curr_node] = val[2][1]
    coords[3, curr_node] = 1.
    curr_node = curr_node + 1

    # face 2
    coords[1, curr_node] = 1.
    coords[2, curr_node] = val[1][1]
    coords[3, curr_node] = val[2][1]
    curr_node = curr_node + 1

    # face 3
    coords[1, curr_node] = val[1][1]
    coords[2, curr_node] = val[2][1]
    coords[3, curr_node] = -1.
    curr_node = curr_node + 1

    # face 4
    coords[1, curr_node] = -1.
    coords[2, curr_node] = val[1][1]
    coords[3, curr_node] = val[2][1]
    curr_node = curr_node + 1

    # face 5
    coords[1, curr_node] = val[1][1]
    coords[2, curr_node] = -1.
    coords[3, curr_node] = val[2][1]
    curr_node = curr_node + 1

    # face 6
    coords[1, curr_node] = val[1][1]
    coords[2, curr_node] = 1.
    coords[3, curr_node] = val[2][1]
    curr_node = curr_node + 1
  end

  # interior nodes
  for val in Iterators.product(edge_coords[2:end - 1], edge_coords[2:end - 1], edge_coords[2:end - 1])
    coords[1, curr_node] = val[1][1]
    coords[2, curr_node] = val[2][1]
    coords[3, curr_node] = val[3][1]
    curr_node = curr_node + 1
  end

  return map(x -> convert_to_vector_coords(e, backend, x...), eachcol(coords))
end

function surface_nodal_coordinates(e::AbstractHex, backend::ArrayBackend)

end

function quadrature_points_and_weights(e::AbstractHex, backend::ArrayBackend)
  ξs, ws = gausslegendre(quadrature_degree(e))
  ξ_return = Matrix{eltype(ξs)}(undef, 3, length(ξs) * length(ξs) * length(ξs))
  w_return = Vector{eltype(ξs)}(undef, length(ξs) * length(ξs) * length(ξs))
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs, ξs))
    ξ_return[1, q] = ξ[1]
    ξ_return[2, q] = ξ[2]
    ξ_return[3, q] = ξ[3]
  end
  for (q, w) in enumerate(Base.Iterators.product(ws, ws, ws))
    w_return[q] = w[1] * w[2] * w[3]
  end
  return map(x -> convert_to(backend, x...), eachcol(ξ_return)), w_return
end

function surface_quadrature_points_and_weights(e::AbstractHex, backend::ArrayBackend)
  ξs, ws = quadrature_points_and_weights(surface_element(e), backend)

  face_1_ξs = map(x -> vcat(x[1], x[2], 1.), ξs)
  face_2_ξs = map(x -> vcat(1., x[1], x[2]), ξs)
  face_3_ξs = map(x -> vcat(x[1], x[2], -1.), ξs)
  face_4_ξs = map(x -> vcat(-1., x[1], x[2]), ξs)
  face_5_ξs = map(x -> vcat(x[1], -1., x[2]), ξs)
  face_6_ξs = map(x -> vcat(x[1], 1., x[2]), ξs)

  return [face_1_ξs, face_2_ξs, face_3_ξs, face_4_ξs, face_5_ξs, face_6_ξs], 
         [ws, ws, ws, ws, ws, ws]
end

# specific implementations

# Constant hex
"""
$(TYPEDEF)
"""
struct Hex0{I, Q} <: AbstractHex{8, I, 0, Q}
end
surface_element(::Hex0{I, Q}) where {I, Q} = Quad0{I, Q}()

# Lagrange implementation
function shape_function_value(e::Hex0{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_vector(e, backend,
    1.
  )
  return Ns
end

function shape_function_gradient(e::Hex0{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_matrix(e, backend,
    0.,
    #
    0.,
    #
    0.
  )
  return Ns
end

function shape_function_hessian(e::Hex0{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_3d_array(e, backend,
    0., 0., 0.,
    #
    0., 0., 0.,
    #
    0., 0., 0.
  )
  return Ns
end

# Linear hex
"""
$(TYPEDEF)
"""
struct Hex8{I, Q} <: AbstractHex{8, I, 1, Q}
end
surface_element(::Hex8{I, Q}) where {I, Q} = Quad4{I, Q}()

# Lagrange implementation
function shape_function_value(e::Hex8{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_vector(e, backend,
    0.125 * (1 - ξ[1]) * (1 - ξ[2]) * (1 - ξ[3]),
    0.125 * (1 + ξ[1]) * (1 - ξ[2]) * (1 - ξ[3]),
    0.125 * (1 + ξ[1]) * (1 + ξ[2]) * (1 - ξ[3]),
    0.125 * (1 - ξ[1]) * (1 + ξ[2]) * (1 - ξ[3]),
    0.125 * (1 - ξ[1]) * (1 - ξ[2]) * (1 + ξ[3]),
    0.125 * (1 + ξ[1]) * (1 - ξ[2]) * (1 + ξ[3]),
    0.125 * (1 + ξ[1]) * (1 + ξ[2]) * (1 + ξ[3]),
    0.125 * (1 - ξ[1]) * (1 + ξ[2]) * (1 + ξ[3]),
  )
  return Ns
end

function shape_function_gradient(e::Hex8{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_matrix(e, backend,
    -0.125 * (1 - ξ[2]) * (1 - ξ[3]),
    0.125 * (1 - ξ[2]) * (1 - ξ[3]),
    0.125 * (1 + ξ[2]) * (1 - ξ[3]),
    -0.125 * (1 + ξ[2]) * (1 - ξ[3]),
    -0.125 * (1 - ξ[2]) * (1 + ξ[3]),
    0.125 * (1 - ξ[2]) * (1 + ξ[3]),
    0.125 * (1 + ξ[2]) * (1 + ξ[3]),
    -0.125 * (1 + ξ[2]) * (1 + ξ[3]),
    #
    -0.125 * (1 - ξ[1]) * (1 - ξ[3]),
    -0.125 * (1 + ξ[1]) * (1 - ξ[3]),
    0.125 * (1 + ξ[1]) * (1 - ξ[3]),
    0.125 * (1 - ξ[1]) * (1 - ξ[3]),
    -0.125 * (1 - ξ[1]) * (1 + ξ[3]),
    -0.125 * (1 + ξ[1]) * (1 + ξ[3]),
    0.125 * (1 + ξ[1]) * (1 + ξ[3]),
    0.125 * (1 - ξ[1]) * (1 + ξ[3]),
    #
    -0.125 * (1 - ξ[1]) * (1 - ξ[2]),
    -0.125 * (1 + ξ[1]) * (1 - ξ[2]),
    -0.125 * (1 + ξ[1]) * (1 + ξ[2]),
    -0.125 * (1 - ξ[1]) * (1 + ξ[2]),
    0.125 * (1 - ξ[1]) * (1 - ξ[2]),
    0.125 * (1 + ξ[1]) * (1 - ξ[2]),
    0.125 * (1 + ξ[1]) * (1 + ξ[2]),
    0.125 * (1 - ξ[1]) * (1 + ξ[2])

  )
  return Ns
end

function shape_function_hessian(e::Hex8{Lagrange}, X, ξ, backend::ArrayBackend)
  Ns = convert_to_3d_array(e, backend,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    #
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    #
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.
  )
  return Ns
end

# general higher order hex
struct Hex{V, I, P, Q} <: AbstractHex{V, I, P, Q}
end
surface_element(::Hex{V, I, P, Q}) where {V, I, P, Q} = Quad{V, I, P, Q}()

function shape_function_value(e::Hex{V, Lagrange}, X, ξ, backend::ArrayBackend) where V
  coords_x = nodal_coordinates(surface_element(surface_element(e)), backend)
  coords_y = nodal_coordinates(surface_element(surface_element(e)), backend)
  coords_z = nodal_coordinates(surface_element(surface_element(e)), backend)
  N_x = shape_function_value(surface_element(surface_element(e)), coords_x, ξ[1], backend)
  N_y = shape_function_value(surface_element(surface_element(e)), coords_y, ξ[2], backend)
  N_z = shape_function_value(surface_element(surface_element(e)), coords_z, ξ[3], backend)
  N_x_rev = reverse(N_x)
  N_y_rev = reverse(N_y)
  N_z_rev = reverse(N_z)
  N = Vector{eltype(coords_x[1])}(undef, num_vertices(e))

  # corner nodes first
  N[1] = N_x[1] * N_y[1] * N_z[1]
  N[2] = N_x[end] * N_y[1] * N_z[1]
  N[3] = N_x[end] * N_y[end] * N_z[1]
  N[4] = N_x[1] * N_y[end] * N_z[1]
  N[5] = N_x[1] * N_y[1] * N_z[end]
  N[6] = N_x[end] * N_y[1] * N_z[end]
  N[7] = N_x[end] * N_y[end] * N_z[end]
  N[8] = N_x[1] * N_y[end] * N_z[end]

  # edge nodes next
  # for n in 2:num_nodes_per_edge(e) - 1
  for n in 2:num_vertices_per_edge(e) - 1
    k = 12 * (n - 2)

    # face 1
    N[8 + k + 1] = N_x[n] * N_y[1] * N_z[1]
    N[8 + k + 2] = N_x[end] * N_y[n] * N_z[1]
    N[8 + k + 3] = N_x_rev[n] * N_y[end] * N_z[1]
    N[8 + k + 4] = N_x[1] * N_y_rev[n] * N_z[1]
    # # face 2
    N[8 + k + 5] = N_x[end] * N_y[1] * N_z[n]
    N[8 + k + 6] = N_x[end] * N_y[n] * N_z[end]
    N[8 + k + 7] = N_x[end] * N_y[end] * N_z_rev[n]
    # face 3
    N[8 + k + 8] = N_x_rev[n] * N_y[1] * N_z[end]
    N[8 + k + 9] = N_x[end] * N_y[n] * N_z[end]
    N[8 + k + 10] = N_x[n] * N_y[end] * N_z[end]
    # N[8 + k + 11] = N_x[1] * N_y_rev[n] * N_z[end]
    # face 4
    N[8 + k + 11] = N_x[1] * N_y[end] * N_z[n]
    N[8 + k + 12] = N_x[1] * N_y[1] * N_z_rev[n]

    # @show 8 + k + 13
    # N[8 + k + 2] = N_x[n] * N_y[n] * N_z[1]
    # #
    # N[8 + k + 3] = N_x[1] * N_y[n] * N_z[1]
    # N[8 + k + 4] = N_x[1] * N_y[end] * N_z[1]
    # face 2
    # N[8 + k + 5] = N_x[n] * N_y[1] * N_z[end]
    # N[8 + k + 6] = N_x[end] * N_y[n] * N_z[end]
    # N[8 + k + 7] = N_x[n] * N_y[end] * N_z[end]
    # N[8 + k + 8] = N_x[1] * N_y[n] * N_z[end]
  end

  curr_node = 8 + 12 * (num_vertices_per_edge(e) - 2)
  for (n, (N_1, N_2)) in enumerate(Iterators.product(N_x[2:end - 1], N_y[2:end - 1]))
    k = curr_node + n
    N[k + 1] = N_1 * N_2 * N_z[end]
    curr_node = curr_node + 1
  end

  for (n, (N_1, N_2)) in enumerate(Iterators.product(N_x[2:end - 1], N_y[2:end - 1]))
    k = curr_node + n
    N[k + 1] = N_x[end] * N_1 * N_2
    curr_node = curr_node + 1
  end

  for (n, (N_1, N_2)) in enumerate(Iterators.product(N_x[2:end - 1], N_y[2:end - 1]))
    k = curr_node + n
    N[k + 1] = N_1 * N_2 * N_z[1]
    curr_node = curr_node + 1
  end

  for (n, (N_1, N_2)) in enumerate(Iterators.product(N_x[2:end - 1], N_y[2:end - 1]))
    k = curr_node + n
    N[k + 1] = N_x[1] * N_1 * N_2
    curr_node = curr_node + 1
  end

  for (n, (N_1, N_2)) in enumerate(Iterators.product(N_x[2:end - 1], N_y[2:end - 1]))
    k = curr_node + n
    N[k + 1] = N_1 * N_y_rev[1] * N_2
    curr_node = curr_node + 1
  end

  for (n, (N_1, N_2)) in enumerate(Iterators.product(N_x[2:end - 1], N_y[2:end - 1]))
    k = curr_node + n
    N[k + 1] = N_1 * N_y[end] * N_2
    curr_node = curr_node + 1
  end

  # now for interior nodes
  m = num_vertices(e) - num_interior_vertices(e) + 1
  for (N_1, N_2, N_3) in Iterators.product(N_x[2:end - 1], N_y[2:end - 1], N_z[2:end - 1])
    N[m] = N_1 * N_2 * N_3
    m = m + 1
  end 
  return convert_to_vector(e, backend, N...)
end

function shape_function_gradient(e::Hex{V, Lagrange}, X, ξ, backend::ArrayBackend) where V
  coords_x = nodal_coordinates(surface_element(surface_element(e)), backend)
  coords_y = nodal_coordinates(surface_element(surface_element(e)), backend)
  coords_z = nodal_coordinates(surface_element(surface_element(e)), backend)
  N_x = shape_function_value(surface_element(surface_element(e)), coords_x, ξ[1], backend)
  N_y = shape_function_value(surface_element(surface_element(e)), coords_y, ξ[2], backend)
  N_z = shape_function_value(surface_element(surface_element(e)), coords_z, ξ[3], backend)
  ∇N_x = shape_function_gradient(surface_element(surface_element(e)), coords_x, ξ[1], backend)
  ∇N_y = shape_function_gradient(surface_element(surface_element(e)), coords_y, ξ[2], backend)
  ∇N_z = shape_function_gradient(surface_element(surface_element(e)), coords_z, ξ[3], backend)

  ∇N = Matrix{eltype(coords_x[1])}(undef, 3, num_vertices(e))

  # corner nodes first
  ∇N[1, 1] = ∇N_x[1] * N_y[1] * N_z[1]
  ∇N[2, 1] = N_x[1] * ∇N_y[1] * N_z[1]
  ∇N[3, 1] = N_x[1] * N_y[1] * ∇N_z[1]
  #
  ∇N[1, 2] = ∇N_x[end] * N_y[1] * N_z[1]
  ∇N[2, 2] = N_x[end] * ∇N_y[1] * N_z[1]
  ∇N[3, 2] = N_x[end] * N_y[1] * ∇N_z[1]
  #
  ∇N[1, 3] = ∇N_x[end] * N_y[end] * N_z[1]
  ∇N[2, 3] = N_x[end] * ∇N_y[end] * N_z[1]
  ∇N[3, 3] = N_x[end] * ∇N_y[end] * ∇N_z[1]
  #
  ∇N[1, 4] = ∇N_x[1] * N_y[end] * N_z[1]
  ∇N[2, 4] = N_x[1] * ∇N_y[end] * N_z[1]
  ∇N[3, 4] = N_x[1] * N_y[end] * ∇N_z[1]
  #
  ∇N[1, 5] = ∇N_x[1] * N_y[1] * N_z[end]
  ∇N[2, 5] = N_x[1] * ∇N_y[1] * N_z[end]
  ∇N[3, 5] = N_x[1] * N_y[1] * ∇N_z[end]
  #
  ∇N[1, 6] = ∇N_x[end] * N_y[1] * N_z[end]
  ∇N[2, 6] = N_x[end] * ∇N_y[1] * N_z[end]
  ∇N[3, 6] = N_x[end] * N_y[1] * ∇N_z[end]
  #
  ∇N[1, 7] = ∇N_x[end] * N_y[end] * N_z[end]
  ∇N[2, 7] = N_x[end] * ∇N_y[end] * N_z[end]
  ∇N[3, 7] = N_x[end] * ∇N_y[end] * ∇N_z[end]
  #
  ∇N[1, 8] = ∇N_x[1] * N_y[end] * N_z[end]
  ∇N[2, 8] = N_x[1] * ∇N_y[end] * N_z[end]
  ∇N[3, 8] = N_x[1] * N_y[end] * ∇N_z[end]
  
  # edge nodes next
  # for n in 2:num_nodes_per_edge(e) - 1
  #   k = 4 * (n - 2)
  #   ∇N[1, 4 + k + 1] = ∇N_x[1] * N_y[n]
  #   ∇N[2, 4 + k + 1] = N_x[1] * ∇N_y[n]
  #   ∇N[1, 4 + k + 2] = ∇N_x[n] * N_y[end]
  #   ∇N[2, 4 + k + 2] = N_x[n] * ∇N_y[end]
  #   ∇N[1, 4 + k + 3] = ∇N_x[end] * N_y[n]
  #   ∇N[2, 4 + k + 3] = N_x[end] * ∇N_y[n]
  #   ∇N[1, 4 + k + 4] = ∇N_x[n] * N_y[1]
  #   ∇N[2, 4 + k + 4] = N_x[n] * ∇N_y[1]
  # end

  # now for interior nodes
  m = num_vertices(e) - num_interior_vertices(e) + 1
  Ns = Iterators.product(N_x[2:end - 1], N_y[2:end - 1])
  ∇Ns = Iterators.product(∇N_x[2:end - 1], ∇N_y[2:end - 1])
  for ((N_1, N_2), (∇N_1, ∇N_2)) in zip(Ns, ∇Ns)
    ∇N[1, m] = ∇N_1 * N_2
    ∇N[2, m] = N_1 * ∇N_2
    m = m + 1
  end 

  return convert_to_matrix(e, backend, ∇N...)
end

function shape_function_hessian(e::Hex{V, Lagrange}, X, ξ, backend::ArrayBackend) where V
  coords_x = nodal_coordinates(surface_element(surface_element(e)), backend)
  coords_y = nodal_coordinates(surface_element(surface_element(e)), backend)
  coords_z = nodal_coordinates(surface_element(surface_element(e)), backend)

  N_x = shape_function_value(surface_element(surface_element(e)), coords_x, ξ[1], backend)
  N_y = shape_function_value(surface_element(surface_element(e)), coords_y, ξ[2], backend)
  N_y = shape_function_value(surface_element(surface_element(e)), coords_z, ξ[3], backend)

  ∇N_x = shape_function_gradient(surface_element(surface_element(e)), coords_x, ξ[1], backend)
  ∇N_y = shape_function_gradient(surface_element(surface_element(e)), coords_y, ξ[2], backend)
  ∇N_z = shape_function_gradient(surface_element(surface_element(e)), coords_z, ξ[3], backend)

  ∇∇N_x = shape_function_hessian(surface_element(surface_element(e)), coords_x, ξ[1], backend)
  ∇∇N_y = shape_function_hessian(surface_element(surface_element(e)), coords_y, ξ[2], backend)
  ∇∇N_z = shape_function_hessian(surface_element(surface_element(e)), coords_z, ξ[3], backend)

  # return N_x * N_y

  ∇∇N = Array{eltype(coords_x[1]), 3}(undef, 3, 3, num_vertices(e))

  # corner nodes first
  # ∇∇N[1, 1, 1] = ∇∇N_x[1] * N_y[1]
  # ∇∇N[1, 2, 1] = ∇N_x[1] * ∇N_y[1]
  # ∇∇N[2, 1, 1] = ∇N_x[1] * ∇N_y[1]
  # ∇∇N[2, 2, 1] = N_x[1] * ∇∇N_y[1]
  # #
  # ∇∇N[1, 1, 2] = ∇∇N_x[end] * N_y[1]
  # ∇∇N[1, 2, 2] = ∇N_x[end] * ∇N_y[1]
  # ∇∇N[2, 1, 2] = ∇N_x[end] * ∇N_y[1]
  # ∇∇N[2, 2, 2] = N_x[end] * ∇∇N_y[1]
  # #
  # ∇∇N[1, 1, 3] = ∇∇N_x[end] * N_y[end]
  # ∇∇N[1, 2, 3] = ∇N_x[end] * ∇N_y[end]
  # ∇∇N[2, 1, 3] = ∇N_x[end] * ∇N_y[end]
  # ∇∇N[2, 2, 3] = N_x[end] * ∇∇N_y[end]
  # #
  # ∇∇N[1, 1, 4] = ∇∇N_x[1] * N_y[end]
  # ∇∇N[1, 2, 4] = ∇N_x[1] * ∇N_y[end]
  # ∇∇N[2, 1, 4] = ∇N_x[1] * ∇N_y[end]
  # ∇∇N[2, 2, 4] = N_x[1] * ∇∇N_y[end]
  
  # edge nodes next
  # for n in 2:num_nodes_per_edge(e) - 1
  #   k = 4 * (n - 2)

  #   ∇∇N[1, 1, 4 + k + 1] = ∇∇N_x[1] * N_y[n]
  #   ∇∇N[1, 2, 4 + k + 1] = ∇N_x[1] * ∇N_y[n]
  #   ∇∇N[2, 1, 4 + k + 1] = ∇N_x[1] * ∇N_y[n]
  #   ∇∇N[2, 2, 4 + k + 1] = N_x[1] * ∇∇N_y[n]
  #   #
  #   ∇∇N[1, 1, 4 + k + 2] = ∇∇N_x[n] * N_y[end]
  #   ∇∇N[1, 2, 4 + k + 2] = ∇N_x[n] * ∇N_y[end]
  #   ∇∇N[2, 1, 4 + k + 2] = ∇N_x[n] * ∇N_y[end]
  #   ∇∇N[2, 2, 4 + k + 2] = N_x[n] * ∇∇N_y[end]
  #   #
  #   ∇∇N[1, 1, 4 + k + 3] = ∇∇N_x[end] * N_y[n]
  #   ∇∇N[1, 2, 4 + k + 3] = ∇N_x[end] * ∇N_y[n]
  #   ∇∇N[2, 1, 4 + k + 3] = ∇N_x[end] * ∇N_y[n]
  #   ∇∇N[2, 2, 4 + k + 3] = N_x[end] * ∇∇N_y[n]
  #   #
  #   ∇∇N[1, 1, 4 + k + 4] = ∇∇N_x[n] * N_y[1]
  #   ∇∇N[1, 2, 4 + k + 4] = ∇N_x[n] * ∇N_y[1]
  #   ∇∇N[2, 1, 4 + k + 4] = ∇N_x[n] * ∇N_y[1]
  #   ∇∇N[2, 2, 4 + k + 4] = N_x[n] * ∇∇N_y[1]
  # end

  # now for interior nodes
  # m = num_nodes(e) - num_interior_nodes(e) + 1
  # Ns = Iterators.product(N_x[2:end - 1], N_y[2:end - 1])
  # ∇Ns = Iterators.product(∇N_x[2:end - 1], ∇N_y[2:end - 1])
  # ∇∇Ns = Iterators.product(∇∇N_x[2:end - 1], ∇∇N_y[2:end - 1])
  # for ((N_1, N_2), (∇N_1, ∇N_2), (∇∇N_1, ∇∇N_2)) in zip(Ns, ∇Ns, ∇∇Ns)
  #   ∇∇N[1, 1, m] = ∇∇N_1 * N_2
  #   ∇∇N[1, 2, m] = ∇N_1 * ∇N_2
  #   ∇∇N[2, 1, m] = ∇N_1 * ∇N_2
  #   ∇∇N[2, 2, m] = N_1 * ∇∇N_2
  #   m = m + 1
  # end 

  return convert_to_3d_array(e, backend, ∇∇N...)
end
