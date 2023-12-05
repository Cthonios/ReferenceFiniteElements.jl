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

# using Hex8(1) as a template
for type_pair in types_to_generate_quadrature(Hex8(Val(1)))
  @eval function setup_hex_quadrature_points!(ξ_return::Vector{$(type_pair[2])}, ξs::T) where T <: AbstractArray
    for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs, ξs))
      ξ_return[q] = $(type_pair[2])(ξ[1], ξ[2], ξ[3])
    end
  end

  """
  """
  @eval function quadrature_points_and_weights(
    e::E, 
    ::Type{$(type_pair[2])}, 
    ::Type{$(type_pair[1])}
  ) where E <: AbstractHex
    ξs, ws = gausslegendre(degree(e))
    n_q_points = length(ws)^3
    ξ_return = Vector{$(type_pair[2])}(undef, n_q_points)
    w_return = Vector{$(type_pair[1])}(undef, n_q_points)
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
