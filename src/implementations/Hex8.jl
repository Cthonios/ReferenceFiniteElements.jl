"""
"""
function element_stencil(::Hex8, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  # first one stupidly causes and allocation
  # @time nodal_coordinates = Ftype[
  #   -1.0  1.0  1.0 -1.0 -1.0  1.0 1.0 -1.0
  #   -1.0 -1.0  1.0  1.0 -1.0 -1.0 1.0  1.0
  #   -1.0 -1.0 -1.0  -1.0 1.0  1.0 1.0  1.0
  # ]
  nodal_coordinates = Matrix{Ftype}(undef, 3, 8)
  nodal_coordinates[1, 1] = -1.0
  nodal_coordinates[2, 1] = -1.0
  nodal_coordinates[3, 1] = -1.0
  #
  nodal_coordinates[1, 2] =  1.0
  nodal_coordinates[2, 2] = -1.0
  nodal_coordinates[3, 2] = -1.0
  #
  nodal_coordinates[1, 3] =  1.0
  nodal_coordinates[2, 3] =  1.0
  nodal_coordinates[3, 3] = -1.0
  #
  nodal_coordinates[1, 4] = -1.0
  nodal_coordinates[2, 4] =  1.0
  nodal_coordinates[3, 4] = -1.0
  #
  nodal_coordinates[1, 5] = -1.0
  nodal_coordinates[2, 5] = -1.0
  nodal_coordinates[3, 5] =  1.0
  #
  nodal_coordinates[1, 6] =  1.0
  nodal_coordinates[2, 6] = -1.0
  nodal_coordinates[3, 6] =  1.0
  #
  nodal_coordinates[1, 7] =  1.0
  nodal_coordinates[2, 7] =  1.0
  nodal_coordinates[3, 7] =  1.0
  #
  nodal_coordinates[1, 8] = -1.0
  nodal_coordinates[2, 8] =  1.0
  nodal_coordinates[3, 8] =  1.0
  #
  # first one stupidly cuases an allocation
  # edge nodes TODO
  # edge_nodes = Matrix{Itype}(undef, 2, 12)
  edge_nodes = Itype[
    1 2 6 5 2 3 7 3 4 8 5 1 
    2 6 5 1 3 7 6 4 8 7 8 4
  ]
  # TODO
  # @time face_nodes = Itype[
  #   1 2 3 1 1 5
  #   2 3 4 5 4 6
  #   6 7 8 8 3 7
  #   5 6 7 4 2 8
  # ]
  face_nodes = Matrix{Itype}(undef, 4, 6)
  face_nodes[1, 1] = 1
  face_nodes[2, 1] = 2
  face_nodes[3, 1] = 6
  face_nodes[4, 1] = 5
  #
  face_nodes[1, 2] = 2
  face_nodes[2, 2] = 3
  face_nodes[3, 2] = 7
  face_nodes[4, 2] = 6
  #
  face_nodes[1, 3] = 3
  face_nodes[2, 3] = 4
  face_nodes[3, 3] = 8
  face_nodes[4, 3] = 7
  #
  face_nodes[1, 4] = 1
  face_nodes[2, 4] = 5
  face_nodes[3, 4] = 8
  face_nodes[4, 4] = 4
  #
  face_nodes[1, 5] = 1
  face_nodes[2, 5] = 4
  face_nodes[3, 5] = 3
  face_nodes[4, 5] = 2
  #
  face_nodes[1, 6] = 5
  face_nodes[2, 6] = 6
  face_nodes[3, 6] = 7
  face_nodes[4, 6] = 8
  #
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, edge_nodes, face_nodes, interior_nodes
end

function shape_function_values(e::Hex8, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SVector, MVector}, A2 <: AbstractArray{<:Number, 1}
}

  N = A1{num_nodes(e), eltype(ξ)}(
    0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]) * (1.0 - ξ[3]),
    0.125 * (1.0 + ξ[1]) * (1.0 - ξ[2]) * (1.0 - ξ[3]),
    0.125 * (1.0 + ξ[1]) * (1.0 + ξ[2]) * (1.0 - ξ[3]),
    0.125 * (1.0 - ξ[1]) * (1.0 + ξ[2]) * (1.0 - ξ[3]),
    0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]) * (1.0 + ξ[3]),
    0.125 * (1.0 + ξ[1]) * (1.0 - ξ[2]) * (1.0 + ξ[3]),
    0.125 * (1.0 + ξ[1]) * (1.0 + ξ[2]) * (1.0 + ξ[3]),
    0.125 * (1.0 - ξ[1]) * (1.0 + ξ[2]) * (1.0 + ξ[3]),
  )
end

"""
"""
function shape_function_gradients(e::Hex8, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SMatrix, MMatrix}, A2 <: AbstractArray{<:Number, 1}
}

  N, D = num_nodes(e), num_dimensions(e)
  ∇N_ξ = A1{N, D, eltype(ξ), N * D}(
    -0.125 * (1.0 - ξ[2]) * (1.0 - ξ[3]),
      0.125 * (1.0 - ξ[2]) * (1.0 - ξ[3]),
      0.125 * (1.0 + ξ[2]) * (1.0 - ξ[3]),
    -0.125 * (1.0 + ξ[2]) * (1.0 - ξ[3]),
    -0.125 * (1.0 - ξ[2]) * (1.0 + ξ[3]),
      0.125 * (1.0 - ξ[2]) * (1.0 + ξ[3]),
      0.125 * (1.0 + ξ[2]) * (1.0 + ξ[3]),
    -0.125 * (1.0 + ξ[2]) * (1.0 + ξ[3]),
    #
    -0.125 * (1.0 - ξ[1]) * (1.0 - ξ[3]),
    -0.125 * (1.0 + ξ[1]) * (1.0 - ξ[3]),
      0.125 * (1.0 + ξ[1]) * (1.0 - ξ[3]),
      0.125 * (1.0 - ξ[1]) * (1.0 - ξ[3]),
    -0.125 * (1.0 - ξ[1]) * (1.0 + ξ[3]),
    -0.125 * (1.0 + ξ[1]) * (1.0 + ξ[3]),
      0.125 * (1.0 + ξ[1]) * (1.0 + ξ[3]),
      0.125 * (1.0 - ξ[1]) * (1.0 + ξ[3]),
    #
    -0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]),
    -0.125 * (1.0 + ξ[1]) * (1.0 - ξ[2]),
    -0.125 * (1.0 + ξ[1]) * (1.0 + ξ[2]),
    -0.125 * (1.0 - ξ[1]) * (1.0 + ξ[2]),
      0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]),
      0.125 * (1.0 + ξ[1]) * (1.0 - ξ[2]),
      0.125 * (1.0 + ξ[1]) * (1.0 + ξ[2]),
      0.125 * (1.0 - ξ[1]) * (1.0 + ξ[2])
  )
end

"""
"""
function shape_function_hessians(e::Hex8, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SArray, MArray}, A2 <: AbstractArray{<:Number, 1}
}

  N, D = num_nodes(e), num_dimensions(e)

  ∇∇N_ξ = A1{Tuple{N, D, D}, eltype(ξ), 3, N * D * D}(
    0.,
    0., 
    0., 
    0., 
    0., 
    0., 
    0., 
    0.,
    #
    0.125 * (1. - ξ[3]),
    -0.125 * (1. - ξ[3]),
    0.125 * (1. - ξ[3]),
    -0.125 * (1. - ξ[3]),
    0.125 * (1. + ξ[3]),
    -0.125 * (1. + ξ[3]),
    0.125 * (1. + ξ[3]),
    -0.125 * (1. + ξ[3]),
    #
    0.125 * (1.0 - ξ[2]),
    -0.125 * (1.0 - ξ[2]),
    -0.125 * (1.0 + ξ[2]),
    0.125 * (1.0 + ξ[2]),
    -0.125 * (1.0 - ξ[2]),
    0.125 * (1.0 - ξ[2]),
    0.125 * (1.0 + ξ[2]),
    -0.125 * (1.0 + ξ[2]),
    #
    0.125 * (1. - ξ[3]),
    -0.125 * (1. - ξ[3]),
    0.125 * (1. - ξ[3]),
    -0.125 * (1. - ξ[3]),
    0.125 * (1. + ξ[3]),
    -0.125 * (1. + ξ[3]),
    0.125 * (1. + ξ[3]),
    -0.125 * (1. + ξ[3]),
    #
    0.,
    0., 
    0., 
    0., 
    0.,
    0., 
    0., 
    0.,
    #
    0.125 * (1.0 - ξ[1]),
    0.125 * (1.0 + ξ[1]),
    -0.125 * (1.0 + ξ[1]),
    -0.125 * (1.0 - ξ[1]),
    -0.125 * (1.0 - ξ[1]),
    -0.125 * (1.0 + ξ[1]),
    0.125 * (1.0 + ξ[1]),
    0.125 * (1.0 - ξ[1]),
    #
    0.125 * (1.0 - ξ[2]),
    -0.125 * (1.0 - ξ[2]),
    -0.125 * (1.0 + ξ[2]),
    0.125 * (1.0 + ξ[2]),
    -0.125 * (1.0 - ξ[2]),
    0.125 * (1.0 - ξ[2]),
    0.125 * (1.0 + ξ[2]),
    -0.125 * (1.0 + ξ[2]),
    #
    0.125 * (1.0 - ξ[1]),
    0.125 * (1.0 + ξ[1]),
    -0.125 * (1.0 + ξ[1]),
    -0.125 * (1.0 - ξ[1]),
    -0.125 * (1.0 - ξ[1]),
    -0.125 * (1.0 + ξ[1]),
    0.125 * (1.0 + ξ[1]),
    0.125 * (1.0 - ξ[1]),
    #
    0.,
    0., 
    0., 
    0., 
    0., 
    0., 
    0., 
    0.
  )
end
