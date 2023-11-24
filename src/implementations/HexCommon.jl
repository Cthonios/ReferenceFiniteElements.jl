"""
"""
abstract type AbstractHex{N, D} <: ReferenceFEType{N, D} end

"""
"""
struct Hex8 <: AbstractHex{8, 3}
  degree::Int
end

# TODO move to a general thing somewhere
types_to_generate = (
  (:(SVector{3, Float32}), :Float32),
  (:(SVector{3, Float64}), :Float64),
  (:(MVector{3, Float32}), :Float32),
  (:(MVector{3, Float64}), :Float64),
  # (:(Vector{Float64}), :Float64)
)
for type_pair in types_to_generate

  @eval function setup_hex_quadrature_points!(ξ_return::Vector{$(type_pair[1])}, ξs::T) where T <: AbstractArray
    for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs, ξs))
      ξ_return[q] = $(type_pair[1])(ξ[1], ξ[2], ξ[3])
    end
  end

  """
  """
  @eval function quadrature_points_and_weights(
    e::E, 
    ::Type{$(type_pair[1])}, 
    ::Type{$(type_pair[2])}
  ) where E <: AbstractHex
    ξs, ws = gausslegendre(degree(e))
    n_q_points = length(ws)^3
    ξ_return = Vector{$(type_pair[1])}(undef, n_q_points)
    w_return = Vector{$(type_pair[2])}(undef, n_q_points)
    setup_hex_quadrature_points!(ξ_return, ξs)
    setup_hex_quadrature_weights!(w_return, ws)
    return ξ_return, w_return
  end
end

"""
"""
function setup_hex_quadrature_weights!(w_return::Vector{Ftype}, ws) where Ftype <: AbstractFloat
  for (q, w) in enumerate(Base.Iterators.product(ws, ws, ws))
    w_return[q] = w[1] * w[2] * w[3]
  end
end
