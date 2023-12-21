# type to aid in making ReferenceFE with StructArrays.jl
"""
Interpolant container for a single quadrature point
"""
struct Interpolants{
  N, D, Ftype, L1, L2, P, Q,
  Vec1 <: AbstractArray{Ftype, 1},
  Vec2 <: AbstractArray{Ftype, 1},
  Mat  <: AbstractArray{Ftype, 2},
  Arr  <: AbstractArray{Ftype, 3}
}
  ξ::Vec1
  w::Ftype
  N::Vec2
  ∇N_ξ::Mat
  ∇∇N_ξ::Arr
end


function setup_interpolants_standard!(Ns, ∇N_ξs, ∇∇N_ξs, e::ReferenceFEType, ξs)
  for (q, ξ) in enumerate(ξs)
    Ns[q]     = shape_function_values(e, SVector, ξ)
    ∇N_ξs[q]  = shape_function_gradients(e, SMatrix, ξ)
    ∇∇N_ξs[q] = shape_function_hessians(e, SArray, ξ)
  end
end

function Interpolants(
  e::ReferenceFEType{N, D, P, Q},
  ::Type{A1}, ::Type{A2}, ::Type{A3}, ::Type{T}
) where {
  N, D, P, Q,
  A1 <: Union{SVector, MVector}, 
  A2 <: Union{SMatrix, MMatrix},
  A3 <: Union{SArray, MArray},
  T <: Number
}

  ξs, ws      = quadrature_points_and_weights(e, A1, T)
  Ns          = Vector{A1{N, T}}(undef, length(ξs))
  ∇N_ξs       = Vector{A2{N, D, T, N * D}}(undef, length(ξs))
  ∇∇N_ξs      = Vector{A3{Tuple{N, D, D}, T, 3, N * D * D}}(undef, length(ξs))

  if typeof(e) <: SimplexTri
    coordinates, _, _ = element_stencil(e, Int64, T)
    setup_interpolants_simplex_tri!(Ns, ∇N_ξs, ∇∇N_ξs, e, ξs, coordinates)
  else
    setup_interpolants_standard!(Ns, ∇N_ξs, ∇∇N_ξs, e, ξs)
  end
  s = StructArray{Interpolants{
    N, D, T, N * D, N * D * D, P, Q,
    eltype(ξs), eltype(Ns), eltype(∇N_ξs), eltype(∇∇N_ξs)
  }}((
    ξ=ξs, w=ws, N=Ns, ∇N_ξ=∇N_ξs, ∇∇N_ξ=∇∇N_ξs
  ))
  return s
end

function Interpolants(
  e::ReferenceFEType{N, D, P, Q},
  ::Type{Vector}, ::Type{Matrix}, ::Type{Array}, ::Type{T}
) where {N, D, P, Q, T  <: Number}

  ξs, ws      = quadrature_points_and_weights(e, Vector, T)
  Ns          = Vector{Vector{T}}(undef, length(ξs))
  ∇N_ξs       = Vector{Matrix{T}}(undef, length(ξs))
  ∇∇N_ξs      = Vector{Array{T, 3}}(undef, length(ξs))
  # setup_interpolants!(Ns, ∇N_ξs, ∇∇N_ξs, e, ξs)
  if typeof(e) <: SimplexTri
    # @assert false
    coordinates, _, _ = element_stencil(e, Int64, T)
    setup_interpolants_simplex_tri!(Ns, ∇N_ξs, ∇∇N_ξs, e, ξs, coordinates)
  else
    setup_interpolants_standard!(Ns, ∇N_ξs, ∇∇N_ξs, e, ξs)
  end
  s = StructArray{Interpolants{
    N, D, T, N * D, N * D * D, P, Q,
    eltype(ξs), eltype(Ns), eltype(∇N_ξs), eltype(∇∇N_ξs)
  }}((
    ξ=ξs, w=ws, N=Ns, ∇N_ξ=∇N_ξs, ∇∇N_ξ=∇∇N_ξs
  ))
  return s
end

function Interpolants{SArray, T}(e::ReferenceFEType) where T
  return Interpolants(e, SVector, SMatrix, SArray, T)
end

function Interpolants{MArray, T}(e::ReferenceFEType) where T
  return Interpolants(e, MVector, MMatrix, MArray, T)
end

function Interpolants{Array, T}(e::ReferenceFEType) where T
  return Interpolants(e, Vector, Matrix, Array, T)
end