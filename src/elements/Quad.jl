struct Quad{I, P, Q} <: AbstractElementType{I, 2, P, Q}
end

num_interior_nodes(e::Quad) = (polynomial_degree(e) - 1) * (polynomial_degree(e) - 1)
num_faces(e::Quad) = 4

function num_nodes(e::Quad) 
  if polynomial_degree(e) < 2
    return 4
  else
    return (polynomial_degree(e) + 1) * (polynomial_degree(e) + 1)
  end
end

num_nodes(e::Quad{I, 0, Q}) where {I, Q} = 4

function num_nodes_per_edge(e::Quad) 
  if polynomial_degree(e) < 2
    return 2
  else
    return polynomial_degree(e) + 1
  end
end

function num_shape_functions(e::Quad{Lagrange, P, Q}) where {P, Q}
  if polynomial_degree(e) == 0
    return 1
  else
    return num_nodes(e)
  end
end

num_nodes_per_face(e::Quad) = 0
surface_element(::Quad{I, P, Q}) where {I, P, Q} = Edge{I, P, Q}()

# TODO
function element_edge_nodes(e::Quad, backend::ArrayBackend)
  edges = Matrix{Int64}(undef, num_nodes_per_edge(e), 4)
  edge_nodes = 1:polynomial_degree(e) + 1

  # corners
  edges[1, 1] = 1
  edges[end, 1] = 2
  edges[1, 2] = 2
  edges[end, 2] = 3
  edges[1, 3] = 3
  edges[end, 3] = 4
  edges[1, 4] = 4
  edges[end, 4] = 1

  # middle edge guys
  curr_node = 5
  if length(edge_nodes) > 2
    for n in 2:length(edge_nodes) - 1
      for m in 1:4
        edges[n, m] = curr_node
        curr_node = curr_node + 1
      end
    end
  end
  return map(x -> convert_to(backend, x...), eachcol(edges))
end

function element_face_nodes(e::Quad, backend::ArrayBackend)
  temp = convert_to(backend, zeros(Int, 0)...)
  face_nodes = [temp]
  return face_nodes
end

# TODO
function element_interior_nodes(e::Quad, backend::ArrayBackend)
  nodes = Vector{Int64}(undef, num_interior_nodes(e))
  for n in axes(nodes, 1)
    nodes[n] = 4 * num_nodes_per_edge(e) - 4 + n
  end
  return nodes
end

function nodal_coordinates(e::Quad, backend::ArrayBackend)
  edge_coords = nodal_coordinates(surface_element(e), backend)
  coords = Matrix{Float64}(undef, 2, num_nodes(e))

  # corner nodes
  coords[1, 1] = -1.
  coords[2, 1] = -1.
  coords[1, 2] =  1.
  coords[2, 2] = -1.
  coords[1, 3] =  1.
  coords[2, 3] =  1.
  coords[1, 4] = -1.
  coords[2, 4] =  1.

  # edge mid point nodes
  curr_node = 5
  @views for val in edge_coords[2:end - 1]
    coords[1, curr_node] = val[1]
    coords[2, curr_node] = -1.
    curr_node = curr_node + 1

    coords[1, curr_node] = 1.
    coords[2, curr_node] = val[1]
    curr_node = curr_node + 1

    coords[1, curr_node] = val[1]
    coords[2, curr_node] = 1.
    curr_node = curr_node + 1

    coords[1, curr_node] = -1.
    coords[2, curr_node] = val[1]
    curr_node = curr_node + 1
  end

  # interior nodes
  for val in Iterators.product(edge_coords[2:end - 1], edge_coords[2:end - 1])
    coords[1, curr_node] = val[1][1]
    coords[2, curr_node] = val[2][1]
    curr_node = curr_node + 1
  end

  return map(x -> convert_to_vector_coords(e, backend, x...), eachcol(coords))
end

function surface_nodal_coordinates(e::Quad{P, Q, I}, backend::ArrayBackend) where {P, Q, I}
  coords = nodal_coordinates(e, backend)
  edges = element_edge_nodes(e, backend)
  # surf_coords = map(x -> coords[x] |> collect, edges)
  surf_coords = mapreduce(x -> coords[x] |> collect, hcat, edges)
  return surf_coords
end

function quadrature_points_and_weights(e::Quad, backend::ArrayBackend)
  ξs, ws = gausslegendre(quadrature_degree(e))
  ξ_return = Matrix{eltype(ξs)}(undef, 2, length(ξs) * length(ξs))
  w_return = Vector{eltype(ξs)}(undef, length(ξs) * length(ξs))
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[1, q] = ξ[1]
    ξ_return[2, q] = ξ[2]
  end
  for (q, w) in enumerate(Base.Iterators.product(ws, ws))
    w_return[q] = w[1] * w[2]
  end
  return map(x -> convert_to(backend, x...), eachcol(ξ_return)), w_return
end

function surface_quadrature_points_and_weights(e::Quad, backend::ArrayBackend)
  ξs, ws = quadrature_points_and_weights(surface_element(e), backend)

  face_1_ξs = map(x -> vcat(x, -1.), ξs)
  face_2_ξs = map(x -> vcat(1., x), ξs)
  face_3_ξs = map(x -> vcat(x, 1.), ξs)
  face_4_ξs = map(x -> vcat(-1., x), ξs)

  return [face_1_ξs, face_2_ξs, face_3_ξs, face_4_ξs], [ws, ws, ws, ws]
end 

function shape_function_value(e::Quad{Lagrange, 0, Q}, X, ξ, backend::ArrayBackend) where Q
  Ns = convert_to_vector(e, backend,
    1.
  )
  return Ns
end

function shape_function_gradient(e::Quad{Lagrange, 0, Q}, X, ξ, backend::ArrayBackend) where Q
  Ns = convert_to_matrix(e, backend,
    0.,
    #
    0.
  )
  return Ns
end

function shape_function_hessian(e::Quad{Lagrange, 0, Q}, X, ξ, backend::ArrayBackend) where Q
  Ns = convert_to_3d_array(e, backend,
    0., 0.,
    #
    0., 0.
  )
  return Ns
end

function shape_function_value(e::Quad{Lagrange, 1, Q}, X, ξ, backend::ArrayBackend{T}) where {Q, T}
  # display(ξ)
  Ns = convert_to_vector(e, backend,
    0.25 * (1 - ξ[1]) * (1 - ξ[2]),
    0.25 * (1 + ξ[1]) * (1 - ξ[2]),
    0.25 * (1 + ξ[1]) * (1 + ξ[2]),
    0.25 * (1 - ξ[1]) * (1 + ξ[2])
  )
  return Ns
end

function shape_function_gradient(e::Quad{Lagrange, 1, Q}, X, ξ, backend::ArrayBackend) where Q
  Ns = convert_to_matrix(e, backend,
    -0.25 * (1 - ξ[2]),
    0.25 * (1 - ξ[2]),
    0.25 * (1 + ξ[2]),
    -0.25 * (1 + ξ[2]),
    #
    -0.25 * (1 - ξ[1]),
    -0.25 * (1 + ξ[1]),
    0.25 * (1 + ξ[1]),
    0.25 * (1 - ξ[1])
  )
  return Ns
end

function shape_function_hessian(e::Quad{Lagrange, 1, Q}, X, ξ, backend::ArrayBackend) where Q
  Ns = convert_to_3d_array(e, backend,
    0., 0.,
    0., 0.,
    0., 0.,
    0., 0.,
    #
    0., 0.,
    0., 0.,
    0., 0.,
    0., 0.
  )
  return Ns
end

for n in 2:25
  @eval function shape_function_value(e::Quad{Lagrange, $n, Q}, X, ξ, backend::ArrayBackend) where Q
    coords_x = nodal_coordinates(surface_element(e), backend)
    coords_y = nodal_coordinates(surface_element(e), backend)
    N_x = shape_function_value(surface_element(e), coords_x, ξ[1], backend)
    N_y = shape_function_value(surface_element(e), coords_y, ξ[2], backend)

    N = Vector{eltype(coords_x[1])}(undef, num_nodes(e))

    # corner nodes first
    N[1] = N_x[1] * N_y[1]
    N[2] = N_x[end] * N_y[1]
    N[3] = N_x[end] * N_y[end]
    N[4] = N_x[1] * N_y[end]

    # edge nodes next
    for n in 2:num_nodes_per_edge(e) - 1
      k = 4 * (n - 2)
      N[4 + k + 1] = N_x[n] * N_y[1]
      N[4 + k + 2] = N_x[end] * N_y[n]
      N[4 + k + 3] = N_x[n] * N_y[end]
      N[4 + k + 4] = N_x[1] * N_y[n]
    end

    # now for interior nodes
    m = num_nodes(e) - num_interior_nodes(e) + 1
    for (N_1, N_2) in Iterators.product(N_x[2:end - 1], N_y[2:end - 1])
      N[m] = N_1 * N_2
      m = m + 1
    end 
    return convert_to_vector(e, backend, N...)
  end

  @eval function shape_function_gradient(e::Quad{Lagrange, $n, Q}, X, ξ, backend::ArrayBackend) where Q
    coords_x = nodal_coordinates(surface_element(e), backend)
    coords_y = nodal_coordinates(surface_element(e), backend)
    N_x = shape_function_value(surface_element(e), coords_x, ξ[1], backend)
    N_y = shape_function_value(surface_element(e), coords_y, ξ[2], backend)
    ∇N_x = shape_function_gradient(surface_element(e), coords_x, ξ[1], backend)
    ∇N_y = shape_function_gradient(surface_element(e), coords_x, ξ[2], backend)

    # return N_x * N_y

    ∇N = Matrix{eltype(coords_x[1])}(undef, 2, num_nodes(e))

    # corner nodes first
    ∇N[1, 1] = ∇N_x[1] * N_y[1]
    ∇N[2, 1] = N_x[1] * ∇N_y[1]
    ∇N[1, 2] = ∇N_x[end] * N_y[1]
    ∇N[2, 2] = N_x[end] * ∇N_y[1]
    ∇N[1, 3] = ∇N_x[end] * N_y[end]
    ∇N[2, 3] = N_x[end] * ∇N_y[end]
    ∇N[1, 4] = ∇N_x[1] * N_y[end]
    ∇N[2, 4] = N_x[1] * ∇N_y[end]
    
    # edge nodes next
    for n in 2:num_nodes_per_edge(e) - 1
      k = 4 * (n - 2)
      ∇N[1, 4 + k + 1] = ∇N_x[1] * N_y[n]
      ∇N[2, 4 + k + 1] = N_x[1] * ∇N_y[n]
      ∇N[1, 4 + k + 2] = ∇N_x[n] * N_y[end]
      ∇N[2, 4 + k + 2] = N_x[n] * ∇N_y[end]
      ∇N[1, 4 + k + 3] = ∇N_x[end] * N_y[n]
      ∇N[2, 4 + k + 3] = N_x[end] * ∇N_y[n]
      ∇N[1, 4 + k + 4] = ∇N_x[n] * N_y[1]
      ∇N[2, 4 + k + 4] = N_x[n] * ∇N_y[1]
    end

    # now for interior nodes
    m = num_nodes(e) - num_interior_nodes(e) + 1
    Ns = Iterators.product(N_x[2:end - 1], N_y[2:end - 1])
    ∇Ns = Iterators.product(∇N_x[2:end - 1], ∇N_y[2:end - 1])
    for ((N_1, N_2), (∇N_1, ∇N_2)) in zip(Ns, ∇Ns)
      ∇N[1, m] = ∇N_1 * N_2
      ∇N[2, m] = N_1 * ∇N_2
      m = m + 1
    end 

    return convert_to_matrix(e, backend, ∇N...)
  end

  @eval function shape_function_hessian(e::Quad{Lagrange, $n, Q}, X, ξ, backend::ArrayBackend) where Q
    coords_x = nodal_coordinates(surface_element(e), backend)
    coords_y = nodal_coordinates(surface_element(e), backend)
    N_x = shape_function_value(surface_element(e), coords_x, ξ[1], backend)
    N_y = shape_function_value(surface_element(e), coords_y, ξ[2], backend)
    ∇N_x = shape_function_gradient(surface_element(e), coords_x, ξ[1], backend)
    ∇N_y = shape_function_gradient(surface_element(e), coords_x, ξ[2], backend)
    ∇∇N_x = shape_function_hessian(surface_element(e), coords_x, ξ[1], backend)
    ∇∇N_y = shape_function_hessian(surface_element(e), coords_x, ξ[2], backend)
    # return N_x * N_y

    ∇∇N = Array{eltype(coords_x[1]), 3}(undef, 2, 2, num_nodes(e))

    # corner nodes first
    ∇∇N[1, 1, 1] = ∇∇N_x[1] * N_y[1]
    ∇∇N[1, 2, 1] = ∇N_x[1] * ∇N_y[1]
    ∇∇N[2, 1, 1] = ∇N_x[1] * ∇N_y[1]
    ∇∇N[2, 2, 1] = N_x[1] * ∇∇N_y[1]
    #
    ∇∇N[1, 1, 2] = ∇∇N_x[end] * N_y[1]
    ∇∇N[1, 2, 2] = ∇N_x[end] * ∇N_y[1]
    ∇∇N[2, 1, 2] = ∇N_x[end] * ∇N_y[1]
    ∇∇N[2, 2, 2] = N_x[end] * ∇∇N_y[1]
    #
    ∇∇N[1, 1, 3] = ∇∇N_x[end] * N_y[end]
    ∇∇N[1, 2, 3] = ∇N_x[end] * ∇N_y[end]
    ∇∇N[2, 1, 3] = ∇N_x[end] * ∇N_y[end]
    ∇∇N[2, 2, 3] = N_x[end] * ∇∇N_y[end]
    #
    ∇∇N[1, 1, 4] = ∇∇N_x[1] * N_y[end]
    ∇∇N[1, 2, 4] = ∇N_x[1] * ∇N_y[end]
    ∇∇N[2, 1, 4] = ∇N_x[1] * ∇N_y[end]
    ∇∇N[2, 2, 4] = N_x[1] * ∇∇N_y[end]
    
    # edge nodes next
    for n in 2:num_nodes_per_edge(e) - 1
      k = 4 * (n - 2)

      ∇∇N[1, 1, 4 + k + 1] = ∇∇N_x[1] * N_y[n]
      ∇∇N[1, 2, 4 + k + 1] = ∇N_x[1] * ∇N_y[n]
      ∇∇N[2, 1, 4 + k + 1] = ∇N_x[1] * ∇N_y[n]
      ∇∇N[2, 2, 4 + k + 1] = N_x[1] * ∇∇N_y[n]
      #
      ∇∇N[1, 1, 4 + k + 2] = ∇∇N_x[n] * N_y[end]
      ∇∇N[1, 2, 4 + k + 2] = ∇N_x[n] * ∇N_y[end]
      ∇∇N[2, 1, 4 + k + 2] = ∇N_x[n] * ∇N_y[end]
      ∇∇N[2, 2, 4 + k + 2] = N_x[n] * ∇∇N_y[end]
      #
      ∇∇N[1, 1, 4 + k + 3] = ∇∇N_x[end] * N_y[n]
      ∇∇N[1, 2, 4 + k + 3] = ∇N_x[end] * ∇N_y[n]
      ∇∇N[2, 1, 4 + k + 3] = ∇N_x[end] * ∇N_y[n]
      ∇∇N[2, 2, 4 + k + 3] = N_x[end] * ∇∇N_y[n]
      #
      ∇∇N[1, 1, 4 + k + 4] = ∇∇N_x[n] * N_y[1]
      ∇∇N[1, 2, 4 + k + 4] = ∇N_x[n] * ∇N_y[1]
      ∇∇N[2, 1, 4 + k + 4] = ∇N_x[n] * ∇N_y[1]
      ∇∇N[2, 2, 4 + k + 4] = N_x[n] * ∇∇N_y[1]
    end

    # now for interior nodes
    m = num_nodes(e) - num_interior_nodes(e) + 1
    Ns = Iterators.product(N_x[2:end - 1], N_y[2:end - 1])
    ∇Ns = Iterators.product(∇N_x[2:end - 1], ∇N_y[2:end - 1])
    ∇∇Ns = Iterators.product(∇∇N_x[2:end - 1], ∇∇N_y[2:end - 1])
    for ((N_1, N_2), (∇N_1, ∇N_2), (∇∇N_1, ∇∇N_2)) in zip(Ns, ∇Ns, ∇∇Ns)
      ∇∇N[1, 1, m] = ∇∇N_1 * N_2
      ∇∇N[1, 2, m] = ∇N_1 * ∇N_2
      ∇∇N[2, 1, m] = ∇N_1 * ∇N_2
      ∇∇N[2, 2, m] = N_1 * ∇∇N_2
      m = m + 1
    end 

    return convert_to_3d_array(e, backend, ∇∇N...)
  end
end