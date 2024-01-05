"""
"""
function element_stencil(::Tet10, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates = Matrix{Ftype}(undef, 3, 10)
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
  nodal_coordinates[1, 5] = 0.5
  nodal_coordinates[2, 5] = 0.0
  nodal_coordinates[3, 5] = 0.0
  #
  nodal_coordinates[1, 6] = 0.5
  nodal_coordinates[2, 6] = 0.5
  nodal_coordinates[3, 6] = 0.0
  #
  nodal_coordinates[1, 7] = 0.0
  nodal_coordinates[2, 7] = 0.5
  nodal_coordinates[3, 7] = 0.0
  #
  nodal_coordinates[1, 8] = 0.0
  nodal_coordinates[2, 8] = 0.0
  nodal_coordinates[3, 8] = 0.5
  #
  nodal_coordinates[1, 9] = 0.5
  nodal_coordinates[2, 9] = 0.0
  nodal_coordinates[3, 9] = 0.5
  #
  nodal_coordinates[1, 10] = 0.0
  nodal_coordinates[2, 10] = 0.5
  nodal_coordinates[3, 10] = 0.5
  #
  # TODO add edge nodes
  edge_nodes = Matrix{Itype}(undef, 2, 6)
  #
  face_nodes = Matrix{Itype}(undef, 6, 4)
  face_nodes[1, 1] = 1
  face_nodes[2, 1] = 5
  face_nodes[3, 1] = 2
  face_nodes[4, 1] = 9
  face_nodes[5, 1] = 4
  face_nodes[6, 1] = 8
  #
  face_nodes[1, 2] = 2
  face_nodes[2, 2] = 6
  face_nodes[3, 2] = 3
  face_nodes[4, 2] = 10
  face_nodes[5, 2] = 4
  face_nodes[6, 2] = 8
  #
  face_nodes[1, 3] = 1
  face_nodes[2, 3] = 8
  face_nodes[3, 3] = 4
  face_nodes[4, 3] = 10
  face_nodes[5, 3] = 3
  face_nodes[6, 3] = 7
  #
  face_nodes[1, 4] = 1
  face_nodes[2, 4] = 7
  face_nodes[3, 4] = 3
  face_nodes[4, 4] = 6
  face_nodes[5, 4] = 2
  face_nodes[6, 4] = 5
  #
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, edge_nodes, face_nodes, interior_nodes
end

"""
"""
function shape_function_values(e::Tet10, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SVector, MVector}, A2 <: AbstractArray{<:Number, 1}
}
  t0 = 1 - ξ[1] - ξ[2] - ξ[3]
  t1 = ξ[1]
  t2 = ξ[2]
  t3 = ξ[3]
  N = A1{num_nodes(e), eltype(ξ)}(
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
end

"""
"""
function shape_function_gradients(e::Tet10, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SMatrix, MMatrix}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  t0 = 1 - ξ[1] - ξ[2] - ξ[3]
  t1 = ξ[1]
  t2 = ξ[2]
  t3 = ξ[3]
  ∇N_ξ = A1{N, D, eltype(ξ), N * D}(
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
end

"""
"""
function shape_function_hessians(e::Tet10, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SArray, MArray}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  ∇∇N_ξ = A1{Tuple{N, D, D}, eltype(ξ), 3, N * D * D}(
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
end
