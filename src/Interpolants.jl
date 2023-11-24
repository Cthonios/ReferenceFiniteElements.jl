# type to aid in making ReferenceFE with StructArrays.jl
"""
Interpolant container for a single quadrature point
"""
struct Interpolants{N, D, Ftype, L1, L2}
  ξ::SVector{D, Ftype}
  w::Ftype
  N::SVector{N, Ftype}
  ∇N_ξ::SMatrix{N, D, Ftype, L1}
  ∇∇N_ξ::SArray{Tuple{N, D, D}, Ftype, 3, L2}
end

"""
Constructor for a StructArray of Interpolants
TODO maybe change the name of this. It might be confusing
"""
function Interpolants(
  e::ReferenceFEType{N, D}, ::Type{Ftype} = Float64
) where {N, D, Ftype}

  ξs_temp, ws = quadrature_points_and_weights(e, SVector{D, Ftype}, Ftype)
  ξs = reinterpret(SVector{D, Ftype}, vec(ξs_temp))
  Ns = Vector{SVector{N, Ftype}}(undef, length(ξs))
  ∇N_ξs = Vector{SMatrix{N, D, Ftype, N * D}}(undef, length(ξs))
  ∇∇N_ξs = Vector{SArray{Tuple{N, D, D}, Ftype, 3, N * D * D}}(undef, length(ξs))
  for (q, ξ) in enumerate(ξs)
    # Ns[q]     = shape_function_values(e, ξ)
    Ns[q]     = shape_function_values(e, SVector, ξ)
    ∇N_ξs[q]  = shape_function_gradients(e, SMatrix, ξ)
    ∇∇N_ξs[q] = shape_function_hessians(e, SArray, ξ)
  end
  s = StructArray{Interpolants{N, D, Ftype, N * D, N * D * D}}((
    ξ=ξs, w=ws, N=Ns, ∇N_ξ=∇N_ξs, ∇∇N_ξ=∇∇N_ξs
  ))
  return s
end
