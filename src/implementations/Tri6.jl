"""
"""
function element_stencil(::Tri6, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates = Ftype[
    0.0 1.0 0.0 0.5 0.5 0.0;
    0.0 0.0 1.0 0.0 0.5 0.5
  ]
  # @time edge_nodes = Itype[
  #   1 2 3
  #   4 5 6
  #   2 3 1
  # ]
  edge_nodes = Matrix{Itype}(undef, 3, 3)
  edge_nodes[1, 1] = 1
  edge_nodes[2, 1] = 4
  edge_nodes[3, 1] = 2
  #
  edge_nodes[1, 2] = 2
  edge_nodes[2, 2] = 5
  edge_nodes[3, 2] = 3
  #
  edge_nodes[1, 1] = 3
  edge_nodes[2, 1] = 6
  edge_nodes[3, 1] = 1
  
  face_nodes = Itype[;;]

  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, edge_nodes, face_nodes, interior_nodes
end

"""
"""
function shape_function_values(e::Tri6, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SVector, MVector}, A2 <: AbstractArray{<:Number, 1}
}
  λ = 1. - ξ[1] - ξ[2]
  N = A1{num_nodes(e), eltype(ξ)}(
    λ * (2. * λ - 1.),
    ξ[1] * (2. * ξ[1] - 1.),
    ξ[2] * (2. * ξ[2] - 1.),
    4. * ξ[1] * λ,
    4. * ξ[1] * ξ[2],
    4. * ξ[2] * λ
  )
end

"""
"""
function shape_function_gradients(e::Tri6, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SMatrix, MMatrix}, A2 <: AbstractArray{<:Number}
}
  λ = 1. - ξ[1] - ξ[2]
  ∇N_ξ = A1{num_nodes(e), num_dimensions(e), eltype(ξ), num_nodes(e) * num_dimensions(e)}(
    -1. * (2. * λ - 1.) - 2. * λ,
     (2. * ξ[1] - 1.) + 2. * ξ[1],
     0.0,
     4. * λ - 4. * ξ[1],
     4. * ξ[2],
    -4. * ξ[2],
    #
    -1. * (2. * λ - 1.) - 2. * λ,
     0.0,
     (2. * ξ[2] - 1.) + 2. * ξ[2],
    -4. * ξ[1],
     4. * ξ[1],
     4. * λ - 4. * ξ[2]
  )
end

"""
"""
function shape_function_hessians(e::Tri6, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SArray, MArray}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  λ = 1. - ξ[1] - ξ[2]
  ∇∇N_ξ = A1{Tuple{N, D, D}, eltype(ξ), 3, N * D * D}(
    4., 4., 0., -8., 0., 0.,
    4., 0., 0., -4., 4., -4.,
    4., 0., 0., -4., 4., -4.,
    4., 0., 4.,  0., 0., -8.
  )
end

