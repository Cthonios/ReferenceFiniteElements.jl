"""
"""
abstract type AbstractQuad{N, D, Q} <: ReferenceFEType{N, D, Q} end

"""
"""
struct Quad4{Q} <: AbstractQuad{4, 2, Q}
  degree::Int64
end

"""
"""
struct Quad9{Q} <: AbstractQuad{9, 2, Q}
  degree::Int64
end

for n in 1:25
  @eval Quad4(::Val{$n}) = Quad4{$(n^2)}($n)
  @eval Quad9(::Val{$n}) = Quad9{$(n^2)}($n)
end

Quad4(q::Int64) = Quad4(Val{q}())
Quad9(q::Int64) = Quad9(Val{q}())

"""
"""
function setup_quad_quadrature_weights!(w_return, ws)
  for (q, w) in enumerate(Base.Iterators.product(ws, ws))
    w_return[q] = w[1] * w[2]
  end
end

function setup_quad_quadrature_points!(ξ_return, ξs)
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[q] = eltype(ξ_return)(ξ[1], ξ[2])
  end
end

"""
"""
function quadrature_points_and_weights(e::E, ::Type{A}, ::Type{T}) where {
  E <: AbstractQuad, A <: Union{SVector, MVector}, T <: Number
}
  D = num_dimensions(e)
  ξs, ws = gausslegendre(degree(e))
  n_q_points = length(ws)^2
  ξ_return = Vector{A{D, T}}(undef, n_q_points)
  w_return = Vector{T}(undef, n_q_points)
  setup_quad_quadrature_points!(ξ_return, ξs)
  setup_quad_quadrature_weights!(w_return, ws)
  return ξ_return, w_return
end
