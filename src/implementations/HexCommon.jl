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

"""
"""
function surface_quadrature_points_and_weights(e::E, ::Type{A}, ::Type{T}) where {
  A <: Union{SVector, MVector}, T <: Number, E <: AbstractHex
}
  D = num_dimensions(e)
  ξs, ws = gausslegendre(degree(e))
  n_s_q_points = length(ws)^2

  ξ_s_return = Matrix{A{D, T}}(undef, 6, n_s_q_points)
  w_s_return = Matrix{T}(undef, 6, n_s_q_points)
  n_s_return = Matrix{A{D, T}}(undef, 6, n_s_q_points)

  setup_hex_surface_quadrature_points!(ξ_s_return, ξs)
  setup_hex_surface_quadrature_weights!(w_s_return, ws)
  setup_hex_surface_normals!(n_s_return, n_s_q_points)
  return ξ_s_return, w_s_return, n_s_return
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

function setup_hex_surface_quadrature_points!(ξ_return, ξs::T) where T <: AbstractArray
  # side 1 front z
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[1, q] = eltype(ξ_return)(ξ[1], ξ[2], 1.0)
  end

  # side 2 right x
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[2, q] = eltype(ξ_return)(1.0, ξ[1], ξ[2])
  end

  # side 3 back z
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[3, q] = eltype(ξ_return)(ξ[1], ξ[2], -1.0)
  end

  # side 4 left x
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[4, q] = eltype(ξ_return)(-1.0, ξ[1], ξ[2])
  end

  # side 5 bottom y
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[5, q] = eltype(ξ_return)(ξ[1], -1.0, ξ[2])
  end

  # side 6 bottom y
  for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
    ξ_return[6, q] = eltype(ξ_return)(ξ[1], -1.0, ξ[2])
  end
end

function setup_hex_surface_quadrature_weights!(w_return, ws)
  for (q, w) in enumerate(Base.Iterators.product(ws, ws))
    for s in 1:6
      w_return[s, q] = w[1] * w[2]
    end
  end
end

function setup_hex_surface_normals!(n_return, n_q_points)
  # side 1 front z
  for q in 1:n_q_points
    n_return[1, q] = eltype(n_return)(0., 0., 1.)
  end

  # side 2 right x
  for q in 1:n_q_points
    n_return[2, q] = eltype(n_return)(1., 0., 0.)
  end

  # side 3 back z
  for q in 1:n_q_points
    n_return[3, q] = eltype(n_return)(0., 0., -1.)
  end

  # side 4 left x
  for q in 1:n_q_points
    n_return[4, q] = eltype(n_return)(-1., 0., 0.)
  end

  # side 5 bottom y
  for q in 1:n_q_points
    n_return[5, q] = eltype(n_return)(0., -1., 0.)
  end

  # side 6 bottom y
  for q in 1:n_q_points
    n_return[6, q] = eltype(n_return)(0., 1., 0.)
  end
end  