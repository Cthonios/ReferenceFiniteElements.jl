struct Quadrature{Rtype <: Real}
  ξs::VecOrMat{Rtype}
  ws::Vector{Rtype}
end
# Quadrature(e::ReferenceFEStencil{I, R, E}) where {I, R, E} = Quadrature(reference_fe_type(e)(), e.degree)
num_dimensions(q::Quadrature{RType}) where RType = size(q.ξs, 1)
num_q_points(q::Quadrature{RType}) where RType = size(q.ξs, 2)
quadrature_points(q::Quadrature{RType}) where RType = q.ξs
quadrature_weights(q::Quadrature{RType}) where RType = q.ws
function Base.show(io::IO, q::Quadrature{RType}) where RType
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
