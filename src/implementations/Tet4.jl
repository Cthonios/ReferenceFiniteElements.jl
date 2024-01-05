"""
"""
function element_stencil(::Tet4, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  # first one stupidly causes an allocation
  # nodal_coordinates = Ftype[
  #   0.0 1.0 0.0 0.0
  #   0.0 0.0 1.0 0.0
  #   0.0 0.0 0.0 1.0
  # ]
  nodal_coordinates = Matrix{Ftype}(undef, 3, 4)
  nodal_coordinates[1, 1] = 0.0
  nodal_coordinates[2, 1] = 0.0
  nodal_coordinates[3, 1] = 0.0
  #
  nodal_coordinates[1, 2] = 1.0
  nodal_coordinates[2, 2] = 0.0
  nodal_coordinates[3, 2] = 0.0
  #
  nodal_coordinates[1, 3] = 0.0
  nodal_coordinates[2, 3] = 1.0
  nodal_coordinates[3, 3] = 0.0
  #
  nodal_coordinates[1, 4] = 0.0
  nodal_coordinates[2, 4] = 0.0
  nodal_coordinates[3, 4] = 1.0
  #
  # TODO add edge nodes
  edge_nodes = Matrix{Itype}(undef, 2, 6)

  # first one stupidly causes an allocation
  # face_nodes = Itype[
  #   1 2 1 1
  #   2 3 4 3
  #   4 4 3 2
  # ]
  face_nodes = Matrix{Itype}(undef, 3, 4)
  face_nodes[1, 1] = 1
  face_nodes[2, 1] = 2
  face_nodes[3, 1] = 4
  #
  face_nodes[1, 2] = 2
  face_nodes[2, 2] = 3
  face_nodes[3, 2] = 4
  #
  face_nodes[1, 3] = 1
  face_nodes[2, 3] = 4
  face_nodes[3, 3] = 3
  #
  face_nodes[1, 4] = 1
  face_nodes[2, 4] = 3
  face_nodes[3, 4] = 2
  #
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, edge_nodes, face_nodes, interior_nodes
end

"""
"""
function shape_function_values(e::Tet4, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SVector, MVector}, A2 <: AbstractArray{<:Number, 1}
}
  N = A1{num_nodes(e), eltype(ξ)}(
    1. - ξ[1] - ξ[2] - ξ[3],
    ξ[1],
    ξ[2],
    ξ[3]
  )
end

"""
"""
function shape_function_gradients(e::Tet4, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SMatrix, MMatrix}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)

  ∇N_ξ = A1{N, D, eltype(ξ), N * D}(
    -1., 1., 0., 0.,
    -1., 0., 1., 0.,
    -1., 0., 0., 1. 
  )
end

"""
"""
function shape_function_hessians(e::Tet4, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SArray, MArray}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  ∇∇N_ξ = zeros(A1{Tuple{N, D, D}, eltype(ξ), 3, N * D * D})
end

