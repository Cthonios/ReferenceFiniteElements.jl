struct Quadrature{Rtype <: Real}
  ξs::VecOrMat{Rtype}
  ws::Vector{Rtype}
end
Quadrature(e::ElementStencil) = Quadrature(e.element_type, e.degree)
num_dimensions(q::Quadrature) = size(q.ξs, 1)
num_q_points(q::Quadrature) = size(q.ξs, 2)
function Base.show(io::IO, q::Quadrature)
  str = "Quadrature:\n"
  str = str * "  Quadrature points = "
  for ξ in eachcol(q.ξs)
    str = str * "$ξ\n                      "
  end
  str = str * "\n"
  str = str * "  Quadrature weights = "
  for w in q.ws
    str = str * "$w\n                       "
  end
  str = str * "\n"
  print(io, str)
end

export Quadrature

export num_dimensions
export num_q_points
