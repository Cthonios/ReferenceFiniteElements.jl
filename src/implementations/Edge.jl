struct Edge <: ReferenceFE
end

function ElementStencil(e::Edge, degree::I, Itype::Type = Integer, Rtype::Type = Float64) where I <: Integer
  degree = degree + 1
  Xs, _ = gausslobatto(degree)
  # to be consistent with optimism
  Xs .= (1. .+ Xs) ./ 2.
  vertex_points = [1, degree]
  face_nodes = Matrix{Integer}(undef, 0, 0)
  interior_points = 2:degree - 1
  return ElementStencil{Itype, Rtype}(
    e, degree - 1, Xs, vertex_points, face_nodes, interior_points
  )
end

function Quadrature(::Edge, degree::I, Rtype::Type = Float64) where I <: Integer
  ξ, w = gausslegendre(degree)
  
  # to make things consistent with optimism
  ξ = (ξ .+ 1.) / 2.
  w = 0.5 * w
  return Quadrature{Rtype}(ξ, w)
end


export Edge
