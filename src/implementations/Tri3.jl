"""
"""
function element_stencil(::Tri3, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates = Ftype[
    0.0 1.0 0.0;
    0.0 0.0 1.0
  ]
  face_nodes = Itype[
    1 2 3
    2 3 1
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, face_nodes, interior_nodes
end

"""
"""
function shape_function_values(e::Tri3, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SVector, MVector}, A2 <: AbstractArray{<:Number, 1}
}
  N = A1{num_nodes(e), eltype(ξ)}(
    1. - ξ[1] - ξ[2],
    ξ[1],
    ξ[2]
  )
end

"""
"""
function shape_function_gradients(e::Tri3, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SMatrix, MMatrix}, A2 <: AbstractArray{<:Number, 1}
}
  ∇N_ξ = A1{num_nodes(e), num_dimensions(e), eltype(ξ), num_nodes(e) * num_dimensions(e)}(
    -1., 1., 0.,
    -1., 0., 1.
  )
end

"""
"""
function shape_function_hessians(e::Tri3, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SArray, MArray}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  ∇∇N_ξ = zeros(A1{Tuple{N, D, D}, eltype(ξ), 3, N * D * D})
end
