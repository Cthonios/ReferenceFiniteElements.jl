struct Edge{I, P, Q} <: AbstractElementType{I, 1, P, Q}
end

num_interior_nodes(e::Edge) = 0
num_faces(e::Edge) = 0
num_nodes(e::Edge) = polynomial_degree(e) + 1
num_nodes_per_edge(e::Edge) = 0
num_nodes_per_face(e::Edge) = 0
num_quadrature_points(::Edge{I, P, Q}) where {I, P, Q} = Q
num_shape_functions(e::Edge) = polynomial_degree(e) == 0 ? 1 : num_nodes(e)
surface_element(::Edge{I, P, Q}) where {I, P, Q} = Vertex{I, P, Q}()

function element_edge_nodes(::Edge, backend::ArrayBackend)
  edge_nodes = Vector{Vector{Int}}(undef, 0)
  return map(x -> convert_to(x, backend), edge_nodes)
end

function element_face_nodes(::Edge, backend::ArrayBackend)
  face_nodes = Vector{Vector{Int}}(undef, 0)
  return map(x -> convert_to(x, backend), face_nodes)
end

function element_interior_nodes(e::Edge, ::ArrayBackend)
  degree = polynomial_degree(e) #+ 1
  return 2:degree - 1 |> collect
end

# TODO need to specialize nodal coordinates for different 
# interpolation rules
function nodal_coordinates(e::Edge, backend::ArrayBackend)
  degree = polynomial_degree(e) + 1
  if polynomial_degree(e) == 0
    nodal_coordinates = [-1., 1.]
  else
    # nodal_coordinates, _ = gausslobatto(degree)
    nodal_coordinates = LinRange(-1., 1., degree) |> collect
  end
  return map(x -> convert_to_vector_coords(e, backend, x), nodal_coordinates) 
end

# TODO need to specialize normals and quadrature for 
# different quadrature rules
function surface_nodal_coordinates(e::Edge, backend::ArrayBackend)
  surface_coordinates = [
    -1.;;
    1.
  ]
  return map(x -> convert_to_vector_coords(e, backend, x...), surface_coordinates)
end

function surface_normals(::Edge, backend::ArrayBackend)
  normals = [[-1.], [1.]]
  return map(x -> convert_to(backend, x...), normals)
end

function quadrature_points_and_weights(e::Edge, backend::ArrayBackend)
  ξs, ws = gausslegendre(quadrature_degree(e))
  return map(x -> convert_to(x, backend), ξs), ws
end

function surface_quadrature_points_and_weights(::Edge, backend::ArrayBackend)
  ξs = map(x -> convert_to(backend, x...), [[-1.], [1.]])
  ws = [[1.], [1.]]
  return ξs, ws
end

# Lagrange
function shape_function_value(e::Edge{Lagrange, 0, Q}, Xs, ξ, backend) where Q
  Ns = convert_to_vector(e, backend,
    1.
  )
  return Ns
end

function shape_function_gradient(e::Edge{Lagrange, 0, Q}, Xs, ξ, backend) where Q
  ∇Ns = convert_to_matrix(e, backend,
    0.
  )
  return ∇Ns
end

function shape_function_hessian(e::Edge{Lagrange, 0, Q}, Xs, ξ, backend) where Q
  ∇∇Ns = convert_to(backend,
    0.
  )
  return ∇∇Ns
end

function shape_function_value(e::Edge{Lagrange, 1, Q}, Xs, ξ, backend) where Q
  Ns = convert_to_vector(e, backend,
    0.5 * (1 - ξ[1]),
    0.5 * (1 + ξ[1])
  )
  return Ns
end

function shape_function_gradient(e::Edge{Lagrange, 1, Q}, Xs, ξ, backend) where Q
  ∇Ns = convert_to_matrix(e, backend,
    -0.5,
    0.5
  )
  return ∇Ns
end

function shape_function_hessian(e::Edge{Lagrange, 1, Q}, Xs, ξ, backend) where Q
  ∇∇Ns = convert_to_3d_array(e, backend,
    0.,
    0.
  )
  return ∇∇Ns
end

function shape_function_value(e::Edge{Lagrange, 2, Q}, Xs, ξ, backend) where Q
  Ns = convert_to_vector(e, backend,
    0.5 * (ξ[1] * ξ[1] - ξ[1]),
    1 - ξ[1] * ξ[1],
    0.5 * (ξ[1] * ξ[1] + ξ[1])
  )
  return Ns
end

function shape_function_gradient(e::Edge{Lagrange, 2, Q}, Xs, ξ, backend) where Q
  Ns = convert_to_matrix(e, backend,
    0.5 * (2 * ξ[1] - 1),
    -2 * ξ[1],
    0.5 * (2 * ξ[1] + 1),
  )
  return Ns
end

function shape_function_hessian(e::Edge{Lagrange, 2, Q}, Xs, ξ, backend) where Q
  Ns = convert_to_3d_array(e, backend,
    1.,
    -2.,
    1.,
  )
  return Ns
end

for n in 3:25
  @eval function shape_function_value(e::Edge{Lagrange, $n, Q}, Xs, ξ, backend::ArrayBackend) where Q
    n_nodes = polynomial_degree(e) + 1
    A = zeros(length(Xs), n_nodes)
    nf = zeros(1, n_nodes)
    for n in 1:polynomial_degree(e) + 1
      p = basis(Legendre, n - 1)
      p = sqrt(2 * (n - 1) + 1) * p
      for m in axes(A, 1)
        A[m, n] = p(Xs[m][1])
      end
      nf[1, n] = p(ξ[1])
    end
    N = A' \ nf'
    return convert_to_vector(e, backend, N[:, 1]...)
  end

  @eval function shape_function_gradient(e::Edge{Lagrange, $n, Q}, Xs, ξ, backend::ArrayBackend) where Q
    n_nodes = polynomial_degree(e) + 1
    A = zeros(length(Xs), n_nodes)
    nf = zeros(1, n_nodes)
    for n in 1:polynomial_degree(e) + 1
      p = basis(Legendre, n - 1)
      p = sqrt(2 * (n - 1) + 1) * p
      for m in axes(A, 1)
        A[m, n] = p(Xs[m][1])
      end
      nf[1, n] = derivative(p)(ξ[1])
    end
    ∇N_ξ = A' \ nf'
    return convert_to_matrix(e, backend, ∇N_ξ[:, 1]...)
  end

  @eval function shape_function_hessian(e::Edge{Lagrange, $n, Q}, Xs, ξ, backend::ArrayBackend) where Q
    n_nodes = polynomial_degree(e) + 1
    A = zeros(length(Xs), n_nodes)
    nf = zeros(1, n_nodes)
    for n in 1:polynomial_degree(e) + 1
      p = basis(Legendre, n - 1)
      p = sqrt(2 * (n - 1) + 1) * p
      for m in axes(A, 1)
        A[m, n] = p(Xs[m][1])
      end
      nf[1, n] = derivative(p, 2)(ξ[1])
    end
    ∇∇N_ξ = A' \ nf'
    return convert_to_3d_array(e, backend, ∇∇N_ξ[:, 1]...)
  end
end
