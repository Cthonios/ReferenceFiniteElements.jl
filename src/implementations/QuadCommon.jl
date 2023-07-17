"""
"""
abstract type AbstractQuad{N, D} <: ReferenceFEType{N, D} end

"""
"""
struct Quad4 <: AbstractQuad{4, 2}
  degree::Int64
end

"""
"""
struct Quad9 <: AbstractQuad{9, 2}
  degree::Int64
end

"""
"""
function quadrature_points_and_weights(e::E, ::Type{Ftype} = Float64) where {E <: AbstractQuad, Ftype <: AbstractFloat}
  ξs, ws = gausslegendre(degree(e))
  n_q_points = length(ws)^2
  ξ_return = Vector{SVector{2, Ftype}}(undef, n_q_points)
  w_return = Vector{Ftype}(undef, n_q_points)
  setup_quad_quadrature_points!(ξ_return, ξs)
  setup_quad_quadrature_weights!(w_return, ws)
  return ξ_return, w_return
end

function setup_quad_quadrature_points!(ξ_return::Vector{SVector{2, Ftype}}, ξs::T) where {Ftype <: AbstractFloat, T <: AbstractArray}
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[q] = SVector{2, Ftype}(ξ[1], ξ[2])
  end
end

function setup_quad_quadrature_weights!(w_return::Vector{Ftype}, ws::T) where {Ftype <: AbstractFloat, T <: AbstractArray}
  for (q, w) in enumerate(Base.Iterators.product(ws, ws))
    w_return[q] = w[1] * w[2]
  end
end
