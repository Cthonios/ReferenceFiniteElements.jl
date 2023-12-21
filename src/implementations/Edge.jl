"""
TODO fix the 1 below and make multiple edge types maybe?
"""
# const Edge = ReferenceFEType{1, 1} where N

abstract type AbstractEdge{N, D, P, Q} <: ReferenceFEType{N, D, P, Q} end

struct Edge2{Q} <: AbstractEdge{2, 1, 1, Q}
  degree::Int64
end

struct Edge3{Q} <: AbstractEdge{3, 1, 2, Q}
  degree::Int64
end

Edge2(::Val{Q}) where Q = Edge2{Q}(Q)
Edge3(::Val{Q}) where Q = Edge3{Q}(Q)

Edge2(q::Int64) = Edge2(Val(q))
Edge3(q::Int64) = Edge3(Val(q))

function element_stencil(e::E, ::Type{Itype}, ::Type{Ftype}) where {
  E <: AbstractEdge, Itype <: Integer, Ftype <: AbstractFloat
}
  # degree = e.degree + 1
  nodal_coordinates, _ = gausslobatto(num_nodes(e))
  # to be consistent with optimism
  # nodal_coordinates .= (1. .+ nodal_coordinates) ./ 2.
  # face_nodes = Matrix{Integer}(undef, 0, 0)
  face_nodes = Itype[1;; num_nodes(e)]
  interior_nodes = 2:num_nodes(e) - 1
  return nodal_coordinates, face_nodes, interior_nodes
end

function quadrature_points_and_weights(e::E, ::Type{A}, ::Type{T}) where {
  A <: Union{SVector, MVector}, T <: Number, E <: AbstractEdge
}

  D = num_dimensions(e)
  ξs, ws = gausslegendre(e.degree)

  ξs_out = Vector{A{D, T}}(undef, num_q_points(e))
  for q in axes(ξs, 1)
    ξs_out[q] = A{D, T}(ξs[q])
  end

  return ξs_out, ws
end

function shape_function_values(e::Edge2, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SVector, MVector}, A2 <: AbstractArray{<:Number, 1}
}
  N = A1{num_nodes(e), eltype(ξ)}(
    0.5 * (1.0 - ξ[1]),
    0.5 * (1.0 + ξ[1])
  )
end

function shape_function_gradients(e::Edge2, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SMatrix, MMatrix}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  ∇N_ξ = A1{N, D, eltype(ξ), N * D}(
    -0.5, 0.5
  )
end

function shape_function_hessians(e::Edge2, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SArray, MArray}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  ∇∇N_ξ = A1{Tuple{N, D, D}, eltype(ξ), 3, N * D * D}(
    0.0, 0.0
  )
end

function shape_function_values(e::Edge3, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SVector, MVector}, A2 <: AbstractArray{<:Number, 1}
}
  N = A1{num_nodes(e), eltype(ξ)}(
    0.25 * (ξ[1]^2 - ξ[1]),
    0.25 * (ξ[1]^2 + ξ[1]),
    1.0 - ξ[1]^2
  )
end

function shape_function_gradients(e::Edge3, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SMatrix, MMatrix}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  ∇N_ξ = A1{N, D, eltype(ξ), N * D}(
    0.50 * (ξ[1] - 1.),
    0.50 * (ξ[1] + 1.),
    -2. * ξ[1]
  )
end

function shape_function_hessians(e::Edge3, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SArray, MArray}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  ∇∇N_ξ = A1{Tuple{N, D, D}, eltype(ξ), 3, N * D * D}(
    0.50, 0.50, -2.0
  )
end