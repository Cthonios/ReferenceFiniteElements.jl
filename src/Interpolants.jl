# type to aid in making ReferenceFE with StructArrays.jl
"""
Interpolant container for a single quadrature point
"""
struct Interpolants{T, N, D, Vec1, Vec2, Vec3, Mat, Arr}
  ξ::Vec1
  w::Vec2 # can be a float or a vector really
  N::Vec3
  ∇N_ξ::Mat
  ∇∇N_ξ::Arr
end

# struct InerpolantsContainer{T, N, D, Vec1, Vec2, Vec3, Mat, Arr}
#   ξ::Vec1
#   w::Vec2
#   N::Vec3
#   ∇N_ξ::Mat
#   ∇∇N_ξ::Arr
# end

function setup_interpolants_standard!(Ns, ∇N_ξs, ∇∇N_ξs, e::ReferenceFEType, ξs)
  for (q, ξ) in enumerate(ξs)
    Ns[q]     = shape_function_values(e, SVector, ξ)
    ∇N_ξs[q]  = shape_function_gradients(e, SMatrix, ξ)
    ∇∇N_ξs[q] = shape_function_hessians(e, SArray, ξ)
  end
end

function Interpolants{StructArray}(
  e::ReferenceFEType{N, D, Q},
  ::Type{A1}, ::Type{A2}, ::Type{A3}, ::Type{T}
) where {
  N, D, Q,
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
    T, N, D,
    eltype(ξs), eltype(Ns), eltype(∇N_ξs), eltype(∇∇N_ξs)
  }}((
    ξ=ξs, w=ws, N=Ns, ∇N_ξ=∇N_ξs, ∇∇N_ξ=∇∇N_ξs
  ))
  return s
end

function Interpolants{StructArray}(
  e::ReferenceFEType{N, D, Q},
  ::Type{Vector}, ::Type{Matrix}, ::Type{Array}, ::Type{T}
) where {N, D, Q, T  <: Number}

  ξs, ws      = quadrature_points_and_weights(e, Vector, T)
  Ns          = Vector{Vector{T}}(undef, length(ξs))
  ∇N_ξs       = Vector{Matrix{T}}(undef, length(ξs))
  ∇∇N_ξs      = Vector{Array{T, 3}}(undef, length(ξs))
  if typeof(e) <: SimplexTri
    coordinates, _, _ = element_stencil(e, Int64, T)
    setup_interpolants_simplex_tri!(Ns, ∇N_ξs, ∇∇N_ξs, e, ξs, coordinates)
  else
    setup_interpolants_standard!(Ns, ∇N_ξs, ∇∇N_ξs, e, ξs)
  end
  s = StructArray{Interpolants{
    T, N, D,
    eltype(ξs), eltype(ws), eltype(Ns), eltype(∇N_ξs), eltype(∇∇N_ξs)
  }}((
    ξ=ξs, w=ws, N=Ns, ∇N_ξ=∇N_ξs, ∇∇N_ξ=∇∇N_ξs
  ))
  return s
end

function Interpolants{Array}(
  e::ReferenceFEType{N, D, Q},
  ::Type{A1}, ::Type{A2}, ::Type{A3}, ::Type{T}
) where {
  N, D, Q,
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

  return Interpolants{
    T, N, D, 
    typeof(ξs), typeof(ws), typeof(Ns), typeof(∇N_ξs), typeof(∇∇N_ξs)
  }(ξs, ws, Ns, ∇N_ξs, ∇∇N_ξs)
end

function Interpolants{Array}(
  e::ReferenceFEType{N, D, Q},
  ::Type{Vector}, ::Type{Matrix}, ::Type{Array}, ::Type{T}
) where {N, D, Q, T  <: Number}

  ξs, ws      = quadrature_points_and_weights(e, Vector, T)
  Ns          = Vector{Vector{T}}(undef, length(ξs))
  ∇N_ξs       = Vector{Matrix{T}}(undef, length(ξs))
  ∇∇N_ξs      = Vector{Array{T, 3}}(undef, length(ξs))
  if typeof(e) <: SimplexTri
    coordinates, _, _ = element_stencil(e, Int64, T)
    setup_interpolants_simplex_tri!(Ns, ∇N_ξs, ∇∇N_ξs, e, ξs, coordinates)
  else
    setup_interpolants_standard!(Ns, ∇N_ξs, ∇∇N_ξs, e, ξs)
  end

  return Interpolants{
    T, N, D, 
    typeof(ξs), typeof(ws), typeof(Ns), typeof(∇N_ξs), typeof(∇∇N_ξs)
  }(ξs, ws, Ns, ∇N_ξs, ∇∇N_ξs)
end

function Interpolants{A, Array, T}(e::ReferenceFEType) where {A, T}
  return Interpolants{A}(e, Vector, Matrix, Array, T)
end

function Interpolants{A, MArray, T}(e::ReferenceFEType) where {A, T}
  return Interpolants{A}(e, MVector, MMatrix, MArray, T)
end

function Interpolants{A, SArray, T}(e::ReferenceFEType) where {A, T}
  return Interpolants{A}(e, SVector, SMatrix, SArray, T)
end

