struct Tri{I, P, Q} <: AbstractElementType{I, 2, P, Q}
end

# num_interior_nodes(e::Tri) = (polynomial_degree(e) - 1) 
num_faces(e::Tri) = 3

function num_nodes(e::Tri)
  if polynomial_degree(e) < 2
    return 3
  else
    return (polynomial_degree(e) + 1) * (polynomial_degree(e) + 2) ÷ 2
  end
end

function num_nodes_per_edge(e::Tri)
  if polynomial_degree(e) < 2
    return 2
  else
    return polynomial_degree(e) + 1
  end
end

function num_shape_functions(e::Tri{Lagrange, P, Q}) where {P, Q}
  if polynomial_degree(e) == 0
    return 1
  else
    return num_nodes(e)
  end
end
num_nodes_per_face(e::Tri) = 0
num_interior_nodes(e::Tri) = num_nodes(e) - num_faces(e) * num_nodes_per_edge(e) + 3
surface_element(::Tri{I, P, Q}) where {I, P, Q} = Edge{I, P, Q}()

function element_edge_nodes(e::Tri, backend)
  edges = Matrix{Int64}(undef, num_nodes_per_edge(e), 3)
  edge_nodes = 1:polynomial_degree(e) + 1

  # corners
  edges[1, 1] = 1
  edges[end, 1] = 2
  edges[1, 2] = 2
  edges[end, 2] = 3
  edges[1, 3] = 3
  edges[end, 3] = 1

  # middle edge guys
  curr_node = 4
  if length(edge_nodes) > 2
    for n in 2:length(edge_nodes) - 1
      for m in 1:3
        edges[n, m] = curr_node
        curr_node = curr_node + 1
      end
    end
  end

  return map(x -> convert_to(backend, x...), eachcol(edges))
end

function element_face_nodes(e::Tri, backend::ArrayBackend)
  temp = convert_to(backend, zeros(Int, 0)...)
  face_nodes = [temp]
  return face_nodes
end

# TODO
function element_interior_nodes(e::Tri, backend::ArrayBackend)
  nodes = Vector{Int64}(undef, num_interior_nodes(e))
  # for n in axes(nodes, 1)
  #   nodes[n] = 4 * num_nodes_per_edge(e) - 4 + n
  # end
  return nodes
end

function nodal_coordinates(e::Tri, backend)
  edge_coords = nodal_coordinates(surface_element(e), backend)
  edge_coords = map(x -> 0.5 * (x .+ 1), edge_coords)
  coords = Matrix{Float64}(undef, 2, num_nodes(e))

  # corner nodes
  # coords[1, 1] = 0.
  # coords[1, 2] = 

  return coords
end

function surface_nodal_coordinates(e::Tri, backend)
  coords = nodal_coordinates(e, backend)
  coords = map(x -> 0.5 * (x .+ 1), coords)
  edges = element_edge_nodes(e, backend)
  surf_coords = map(x -> coords[x] |> collect, edges)
  return surf_coords
end

