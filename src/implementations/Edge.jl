"""
TODO fix the 1 below and make multiple edge types maybe?
"""
const Edge = ReferenceFEType{1, 1} where N

function element_stencil(::Edge, degree::I, ::Type{Itype}, ::Type{Ftype}) where {I <: Integer, Itype <: Integer, Ftype <: AbstractFloat}
  degree = degree + 1
  nodal_coordinates, _ = gausslobatto(degree)
  # to be consistent with optimism
  nodal_coordinates .= (1. .+ nodal_coordinates) ./ 2.
  edge_nodes = Itype[;;]
  face_nodes = Itype[;;]
  interior_nodes = 2:degree - 1
  return nodal_coordinates, edge_nodes, face_nodes, interior_nodes
end

function quadrature_points_and_weights(::Edge, degree::I, ::Type{Ftype} = Float64) where {I <: Integer, Ftype <: AbstractFloat}
  ξs, ws = gausslegendre(degree)
  
  # to make things consistent with optimism
  ξs .= (ξs .+ 1.) ./ 2.
  ws .= 0.5 .* ws
  return ξs, ws
end
