struct Vertex{I, P, Q} <: AbstractVertex{I, P, Q}
end

num_interior_nodes(e::AbstractVertex) = 0
num_nodes(e::AbstractVertex) = 1
num_nodes_per_Vertex(e::AbstractVertex) = 0
num_nodes_per_face(e::AbstractVertex) = 0
num_quadrature_points(e::AbstractVertex) = 0
surface_element(::AbstractVertex) = nothing

function element_edge_nodes(::Vertex, backend::ArrayBackend)
  Vertex_nodes = Vector{Vector{Int}}(undef, 0)
  return map(x -> convert_to(x, backend), Vertex_nodes)
end

function element_face_nodes(::Vertex, backend::ArrayBackend)
  face_nodes = Vector{Vector{Int}}(undef, 0)
  return map(x -> convert_to(x, backend), face_nodes)
end

function element_interior_nodes(::Vertex, backend::ArrayBackend)
  interior_nodes = Vector{Int}(undef, 0)
  return interior_nodes
end

function nodal_coordinates(::Vertex, backend::ArrayBackend)
  nodal_coordinates = [[0.]]
  return map(x -> convert_to(x, backend), nodal_coordinates) 
end

function surface_nodal_coordinates(::Vertex, backend::ArrayBackend)
  surface_nodal_coordinates = [[0.;;]]
  return map(x -> convert_to(backend, x...), surface_nodal_coordinates) 
end

function surface_normals(::Vertex, backend)
  return map(x -> convert_to(backend, x...), [[;;]])
end

function quadrature_points_and_weights(::Vertex, backend::ArrayBackend)
  ξs = [[0.]]
  ws = [1.]
  ξs = map(x -> convert_to(backend, x...), ξs)
  ws = convert_to(backend, ws)
  return ξs, ws
end

function surface_quadrature_points_and_weights(::Vertex, backend::ArrayBackend)
  surface_ξs = [[;;]]
  surface_ws = [[]]
  surface_ξs = map(x -> convert_to(backend, x...), surface_ξs)
  surface_ws = convert_to(backend, surface_ws)
  return surface_ξs, surface_ws
end

function shape_function_value(::Vertex, X, ξ, backend)
  return convert_to(backend, 1.)
end

function shape_function_gradient(::Vertex, X, ξ, backend)
  return convert_to(backend, 0.)
end

function shape_function_hessian(::Vertex, X, ξ, backend)
  return convert_to(backend, 0.)
end