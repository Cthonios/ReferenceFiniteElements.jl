# # abstract methods
# num_faces(e::AbstractTet) = 4

# function num_nodes(e::AbstractTet)
#   if polynomial_degree(e) < 2
#     return 3
#   else
#     return (polynomial_degree(e) + 1) * (polynomial_degree(e) + 2) * (polynomial_degree(e) + 2) ÷ 6
#   end
# end

# function num_nodes_per_edge(e::AbstractTet)
#   if polynomial_degree(e) < 2
#     return 2
#   else
#     return polynomial_degree(e) + 1
#   end
# end

# function num_shape_functions(e::AbstractTet{Lagrange, P, Q}) where {P, Q}
#   if polynomial_degree(e) == 0
#     return 1
#   else
#     return num_nodes(e)
#   end
# end
# num_nodes_per_face(e::AbstractTet) = num_nodes(surface_element(e))
# num_interior_nodes(e::AbstractTet{I, 0, Q}) where {I, Q} = 0
# num_interior_nodes(e::AbstractTet) = num_nodes(e) - num_faces(e) * num_nodes_per_edge(e) + 4
# surface_element(::AbstractTet{I, P, Q}) where {I, P, Q} = Tri{I, P, Q}()

# # TODO
# function element_edge_nodes(e::AbstractTri, backend)
#   edges = Matrix{Int64}(undef, num_nodes_per_edge(e), 3)
#   edge_nodes = 1:polynomial_degree(e) + 1

#   # corners
#   # edges[1, 1] = 1
#   # edges[end, 1] = 2
#   # edges[1, 2] = 2
#   # edges[end, 2] = 3
#   # edges[1, 3] = 3
#   # edges[end, 3] = 1

#   # # middle edge guys
#   # curr_node = 4
#   # if length(edge_nodes) > 2
#   #   for n in 2:length(edge_nodes) - 1
#   #     for m in 1:3
#   #       edges[n, m] = curr_node
#   #       curr_node = curr_node + 1
#   #     end
#   #   end
#   # end

#   return map(x -> convert_to(backend, x...), eachcol(edges))
# end

# function element_face_nodes(e::AbstractTri, backend::ArrayBackend)
#   temp = convert_to(backend, zeros(Int, 0)...)
#   face_nodes = [temp]
#   return face_nodes
# end

# # TODO
# function element_interior_nodes(e::AbstractTri, backend::ArrayBackend)
#   nodes = Vector{Int64}(undef, num_interior_nodes(e))
#   # for n in axes(nodes, 1)
#   #   nodes[n] = 4 * num_nodes_per_edge(e) - 4 + n
#   # end
#   return nodes
# end

# function nodal_coordinates(e::AbstractTet, backend::ArrayBackend)
#   edge_coords = nodal_coordinates(surface_element(e), backend)
#   edge_coords = map(x -> 0.5 * (x .+ 1), edge_coords)
#   coords = Matrix{Float64}(undef, 3, num_nodes(e))

#   # corner nodes
#   coords[1, 1] = 0.
#   coords[2, 1] = 0.
#   coords[3, 1] = 0.
#   #
#   coords[1, 2] = 1.
#   coords[2, 2] = 0.
#   coords[3, 2] = 0.
#   #
#   coords[1, 3] = 0.
#   coords[2, 3] = 1.
#   coords[3, 3] = 0.
#   #
#   coords[1, 4] = 0.
#   coords[2, 4] = 0.
#   coords[3, 4] = 1.
  
#   # edge mid points

#   # cell interior points
#   return map(x -> convert_to_vector_coords(e, backend, x...), eachcol(coords))
# end

# function surface_nodal_coordinates(e::AbstractTet, backend::ArrayBackend)
#   coords = surface_nodal_coordinates(surface_element(e), backend)
#   edges = element_edge_nodes(e, backend)
  
# end

# function quadrature_points_and_weights(e::AbstractTet, backend::ArrayBackend)
#   if quadrature_degree(e) == 1
#     ξs = Matrix{Float64}(undef, 3, 1)
#     ξs[:, 1] = [1. / 4., 1. / 4., 1. / 4.]
#     ws = [1. / 6.]
#   elseif quadrature_degree(e) == 2
#     ξs = Matrix{Float64}(undef, 3, 5)
#     ξs[:, 1] = [1. / 4., 1. / 4., 1. / 4.]
#     ξs[:, 2] = [1. / 6., 1. / 6., 1. / 6.]
#     ξs[:, 3] = [1. / 6., 1. / 6., 1. / 2.]
#     ξs[:, 4] = [1. / 6., 1. / 2., 1. / 6.]
#     ξs[:, 5] = [1. / 2., 1. / 6., 1. / 6.]

#     #
#     ws = [
#       -2. / 15.
#        3. / 40.
#        3. / 40.
#        3. / 40.
#        3. / 40.
#     ]
#   end
#   return map(x -> convert_to(backend, x...), eachcol(ξs)), ws
# end

# function surface_quadrature_points_and_weights(e::AbstractTet, backend::ArrayBackend)
#   ξs, ws = quadrature_points_and_weights(surface_element(e), backend)

#   face_1_ξs = map(x -> vcat(x, 0., 0.), ξs)
#   face_2_ξs = map(x -> vcat(x, 1. .- x), ξs)
#   face_3_ξs = map(x -> vcat(0., x), ξs)

#   return [face_1_ξs, face_2_ξs, face_3_ξs], [ws, ws, ws]
# end 