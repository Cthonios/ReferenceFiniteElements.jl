"""
"""
const Quad4 = ReferenceFEType{4, 2}
"""
"""
const Quad9 = ReferenceFEType{9, 2}
"""
"""
const QuadUnion = Union{Quad4, Quad9}

"""
"""
function quadrature_points_and_weights(::E, degree::I, ::Type{Ftype} = Float64) where {E <: QuadUnion, I <: Integer, Ftype <: AbstractFloat}
  ξs, ws = gausslegendre(degree)
  n_q_points = length(ws)^2
  ξ_return = zeros(Ftype, 2, n_q_points)
  w_return = zeros(Ftype, n_q_points)
  setup_quad_quadrature_points!(ξ_return, ξs)
  setup_quad_quadrature_weights!(w_return, ws)
  return ξ_return, w_return
end

function setup_quad_quadrature_points!(ξ_return::Matrix{Ftype}, ξs) where Ftype <: AbstractFloat
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[1, q] = ξ[1]
    ξ_return[2, q] = ξ[2]
  end
end

function setup_quad_quadrature_weights!(w_return::Vector{Ftype}, ws) where Ftype <: AbstractFloat
  for (q, w) in enumerate(Base.Iterators.product(ws, ws))
    w_return[q] = w[1] * w[2]
  end
end
