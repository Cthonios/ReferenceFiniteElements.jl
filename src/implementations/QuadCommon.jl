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

"""
"""
function surface_quadrature_points_and_weights(e::E, ::Type{A}, ::Type{T}) where {
  A <: Union{SVector, MVector}, T <: Number, E <: AbstractQuad
}
  D = num_dimensions(e)
  ξs, ws = gausslegendre(degree(e))
  n_s_q_points = length(ws)

  ξ_s_return = Matrix{A{D, T}}(undef, 4, n_s_q_points)
  w_s_return = Matrix{T}(undef, 4, n_s_q_points)
  n_s_return = Matrix{A{D, T}}(undef, 4, n_s_q_points)

  setup_quad_surface_quadrature_points!(ξ_s_return, ξs)
  setup_quad_surface_quadrature_weights!(w_s_return, ws)
  setup_quad_surface_normals!(n_s_return, n_s_q_points)
  return ξ_s_return, w_s_return, n_s_return
end

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

function setup_quad_surface_quadrature_points!(ξ_return, ξs::T) where T <: AbstractArray
  # side 1 negative y
  for (q, ξ) in enumerate(ξs)
    ξ_return[1, q] = eltype(ξ_return)(ξ, -1.0)
  end

  # side 2 positive x
  for (q, ξ) in enumerate(ξs)
    ξ_return[2, q] = eltype(ξ_return)(1.0, ξ)
  end

  # side 3 positive y
  for (q, ξ) in enumerate(ξs)
    ξ_return[3, q] = eltype(ξ_return)(ξ, 1.0)
  end

  # side 4 negative x
  for (q, ξ) in enumerate(ξs)
    ξ_return[4, q] = eltype(ξ_return)(-1.0, ξ)
  end
end

function setup_quad_surface_quadrature_weights!(w_return, ws)
  for (q, w) in enumerate(ws)
    for s in 1:4
      w_return[s, q] = w[1]
    end
  end
end

function setup_quad_surface_normals!(n_return, n_q_points)
  # side 1 negative y
  for q in 1:n_q_points
    n_return[1, q] = eltype(n_return)(0., -1.)
  end

  # side 2 positive x
  for q in 1:n_q_points
    n_return[2, q] = eltype(n_return)(1., 0.)
  end

  # side 3 positive y 
  for q in 1:n_q_points
    n_return[3, q] = eltype(n_return)(0., 1.)
  end

  # side 4 negative x
  for q in 1:n_q_points
    n_return[4, q] = eltype(n_return)(-1., 0.)
  end
end  