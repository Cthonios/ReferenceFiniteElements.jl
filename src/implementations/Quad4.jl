"""
"""
function element_stencil(::Quad4, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates = Ftype[
    -1.0  1.0 1.0 -1.0;
    -1.0 -1.0 1.0  1.0
  ]
  edge_nodes = Itype[
    1 2 3 4
    2 3 4 1
  ]
  face_nodes = Itype[
    ;;
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, edge_nodes, face_nodes, interior_nodes
end

"""
"""
function shape_function_values(e::Quad4, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SVector, MVector}, A2 <: AbstractArray{<:Number, 1}
}
  N = A1{num_nodes(e), eltype(ξ)}(
    0.25 * (1.0 - ξ[1]) * (1.0 - ξ[2]),
    0.25 * (1.0 + ξ[1]) * (1.0 - ξ[2]),
    0.25 * (1.0 + ξ[1]) * (1.0 + ξ[2]),
    0.25 * (1.0 - ξ[1]) * (1.0 + ξ[2])
  )
end

"""
"""
function shape_function_gradients(e::Quad4, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SMatrix, MMatrix}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  ∇N_ξ = A1{N, D, eltype(ξ), N * D}(
    -0.25 * (1.0 - ξ[2]),
      0.25 * (1.0 - ξ[2]),
      0.25 * (1.0 + ξ[2]),
    -0.25 * (1.0 + ξ[2]),
    #
    -0.25 * (1.0 - ξ[1]),
    -0.25 * (1.0 + ξ[1]),
      0.25 * (1.0 + ξ[1]),
      0.25 * (1.0 - ξ[1])
  )
end

"""
"""
function shape_function_hessians(e::Quad4, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SArray, MArray}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  ∇∇N_ξ = A1{Tuple{N, D, D}, eltype(ξ), 3, N * D * D}(
    0.0, 0.0, 0.0, 0.0,
    0.25, -0.25, 0.25, -0.25,
    #
    0.25, -0.25, 0.25, -0.25,
    0.0, 0.0, 0.0, 0.0
  )
end
