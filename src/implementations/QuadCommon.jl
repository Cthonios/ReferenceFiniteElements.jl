"""
"""
abstract type AbstractQuad <: AbstractReferenceFE end

"""
"""
function Quadrature(::E, degree::I, Rtype::Type = Float64) where {I <: Integer, E <: AbstractQuad}
  ξs, ws = gausslegendre(degree)
  n_q_points = length(ws)^2
  ξ_return = zeros(Rtype, 2, n_q_points)
  w_return = zeros(Rtype, n_q_points)
  setup_quad_quadrature_points!(ξ_return, ξs)
  setup_quad_quadrature_weights!(w_return, ws)
  return Quadrature(ξ_return, w_return)
end

function setup_quad_quadrature_points!(ξ_return::Matrix{Rtype}, ξs) where Rtype <: Real
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[1, q] = ξ[1]
    ξ_return[2, q] = ξ[2]
  end
end

function setup_quad_quadrature_weights!(w_return::Vector{Rtype}, ws) where Rtype <: Real
  for (q, w) in enumerate(Base.Iterators.product(ws, ws))
    w_return[q] = w[1] * w[2]
  end
end
