"""
"""
const Hex8 = ReferenceFEType{8, 3}
const HexUnion = Union{Hex8}

"""
"""
function Quadrature(::E, degree::I, Rtype::Type = Float64) where {I <: Integer, E <: HexUnion}
  ξs, ws = gausslegendre(degree)
  n_q_points = length(ws)^3
  ξ_return = zeros(Rtype, 3, n_q_points)
  w_return = zeros(Rtype, n_q_points)
  setup_hex_quadrature_points!(ξ_return, ξs)
  setup_hex_quadrature_weights!(w_return, ws)
  return Quadrature(ξ_return, w_return)
end

function setup_hex_quadrature_points!(ξ_return::Matrix{Rtype}, ξs) where Rtype <: Real
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs, ξs))
    ξ_return[1, q] = ξ[1]
    ξ_return[2, q] = ξ[2]
    ξ_return[3, q] = ξ[3]
  end
end

function setup_hex_quadrature_weights!(w_return::Vector{Rtype}, ws) where Rtype <: Real
  for (q, w) in enumerate(Base.Iterators.product(ws, ws, ws))
    w_return[q] = w[1] * w[2] * w[3]
  end
end
