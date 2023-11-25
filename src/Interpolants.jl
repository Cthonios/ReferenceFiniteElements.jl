# type to aid in making ReferenceFE with StructArrays.jl
"""
Interpolant container for a single quadrature point
"""
struct Interpolants{
  N, D, Ftype, L1, L2,
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

types_to_generate = (
  (
    :(Float32),
    :(SVector),
    :(SMatrix),
    :(SArray)
  ),
  (
    :(Float32),
    :(MVector),
    :(MMatrix),
    :(MArray)
  ),
  (
    :(Float64),
    :(SVector),
    :(SMatrix),
    :(SArray)
  ),
  (
    :(Float64),
    :(MVector),
    :(MMatrix),
    :(MArray)
  )
)

for type in types_to_generate
  @eval function Interpolants(
    e::ReferenceFEType{N, D},
    ::Type{$(type[4])},
    ::Type{$(type[1])}
  ) where {N, D}

    ξs_temp, ws = quadrature_points_and_weights(e, $(type[2]){D, $(type[1])}, $(type[1]))
    ξs          = reinterpret($(type[2]){D, $(type[1])}, vec(ξs_temp))
    Ns          = Vector{$(type[2]){N, $(type[1])}}(undef, length(ξs))
    ∇N_ξs       = Vector{$(type[3]){N, D, $(type[1]), N * D}}(undef, length(ξs))
    ∇∇N_ξs      = Vector{$(type[4]){Tuple{N, D, D}, $(type[1]), 3, N * D * D}}(undef, length(ξs))
    for (q, ξ) in enumerate(ξs)
      Ns[q]     = shape_function_values(e, SVector, ξ)
      ∇N_ξs[q]  = shape_function_gradients(e, SMatrix, ξ)
      ∇∇N_ξs[q] = shape_function_hessians(e, SArray, ξ)
    end
    s = StructArray{Interpolants{
      N, D, $(type[1]), N * D, N * D * D,
      eltype(ξs), eltype(Ns), eltype(∇N_ξs), eltype(∇∇N_ξs)
    }}((
      ξ=ξs, w=ws, N=Ns, ∇N_ξ=∇N_ξs, ∇∇N_ξ=∇∇N_ξs
    ))
    return s
  end
end
