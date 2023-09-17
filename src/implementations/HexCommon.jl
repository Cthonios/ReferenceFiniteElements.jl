"""
"""
abstract type AbstractHex{N, D} <: ReferenceFEType{N, D} end

"""
"""
struct Hex8 <: AbstractHex{8, 3}
  degree::Int
end

"""
"""
function quadrature_points_and_weights(e::E, ::Type{Ftype} = Float64) where {E <: AbstractHex, Ftype <: AbstractFloat}
  ξs, ws = gausslegendre(degree(e))
  n_q_points = length(ws)^3
  ξ_return = Vector{SVector{3, Ftype}}(undef, n_q_points)
  w_return = Vector{Ftype}(undef, n_q_points)
  setup_hex_quadrature_points!(ξ_return, ξs)
  setup_hex_quadrature_weights!(w_return, ws)
  return ξ_return, w_return
end

# TODO could maybe do reinterpret on the below
"""
"""
function setup_hex_quadrature_points!(ξ_return::Vector{SVector{3, Ftype}}, ξs::T) where {Ftype <: AbstractFloat, T <: AbstractArray}
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs, ξs))
    ξ_return[q] = SVector{3, Ftype}(ξ[1], ξ[2], ξ[3])
  end
end

"""
"""
function setup_hex_quadrature_weights!(w_return::Vector{Ftype}, ws) where Ftype <: AbstractFloat
  for (q, w) in enumerate(Base.Iterators.product(ws, ws, ws))
    w_return[q] = w[1] * w[2] * w[3]
  end
end
