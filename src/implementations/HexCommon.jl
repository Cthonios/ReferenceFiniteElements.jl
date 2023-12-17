"""
"""
abstract type AbstractHex{N, D, Q} <: ReferenceFEType{N, D, Q} end

"""
"""
struct Hex8{Q} <: AbstractHex{8, 3, Q}
  degree::Int
end

for n in 1:25
  @eval Hex8(::Val{$n}) = Hex8{$(n^3)}($n)
end 
Hex8(q::Int64) = Hex8(Val{q}())

"""
"""
function quadrature_points_and_weights(e::E, ::Type{A}, ::Type{T}) where {
  A <: Union{SVector, MVector}, T <: Number, E <: AbstractHex
}
  D = num_dimensions(e)
  ξs, ws = gausslegendre(degree(e))
  n_q_points = length(ws)^3
  ξ_return = Vector{A{D, T}}(undef, n_q_points)
  w_return = Vector{T}(undef, n_q_points)
  setup_hex_quadrature_points!(ξ_return, ξs)
  setup_hex_quadrature_weights!(w_return, ws)
  return ξ_return, w_return
end

function setup_hex_quadrature_points!(ξ_return, ξs::T) where T <: AbstractArray
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs, ξs))
    ξ_return[q] = eltype(ξ_return)(ξ[1], ξ[2], ξ[3])
  end
end

function setup_hex_quadrature_weights!(w_return, ws)
  for (q, w) in enumerate(Base.Iterators.product(ws, ws, ws))
    w_return[q] = w[1] * w[2] * w[3]
  end
end
