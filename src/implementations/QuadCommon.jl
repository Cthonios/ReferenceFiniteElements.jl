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


types_to_generate = (
  (:(SVector{2, Float32}), :Float32),
  (:(SVector{2, Float64}), :Float64),
  (:(MVector{2, Float32}), :Float32),
  (:(MVector{2, Float64}), :Float64),
  # (:(Vector{Float64}), :Float64)
)

for type_pair in types_to_generate

  @eval function setup_quad_quadrature_points!(ξ_return::Vector{$(type_pair[1])}, ξs::T) where T <: AbstractArray
    for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
      ξ_return[q] = $(type_pair[1])(ξ[1], ξ[2])
    end
  end

  """
  """
  @eval function quadrature_points_and_weights(
    e::E, 
    ::Type{$(type_pair[1])}, 
    ::Type{$(type_pair[2])}
  ) where E <: AbstractQuad
    ξs, ws = gausslegendre(degree(e))
    n_q_points = length(ws)^2
    ξ_return = Vector{$(type_pair[1])}(undef, n_q_points)
    w_return = Vector{$(type_pair[2])}(undef, n_q_points)
    setup_quad_quadrature_points!(ξ_return, ξs)
    setup_quad_quadrature_weights!(w_return, ws)
    return ξ_return, w_return
  end
end

"""
"""
function setup_quad_quadrature_weights!(w_return::Vector{Ftype}, ws::T) where {Ftype <: AbstractFloat, T <: AbstractArray}
  for (q, w) in enumerate(Base.Iterators.product(ws, ws))
    w_return[q] = w[1] * w[2]
  end
end
