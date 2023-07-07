struct ReferenceFEStencil{IType, RType, RefFE <: ReferenceFEType}
  degree::Integer
  coordinates::VecOrMat{RType}
  vertex_nodes::Vector{IType}
  face_nodes::Matrix{IType}
  interior_nodes::Vector{IType}
end
reference_fe_coordinates(e::ReferenceFEStencil) = e.coordinates
reference_fe_degree(e::ReferenceFEStencil) = e.degree
reference_fe_face_nodes(e::ReferenceFEStencil) = e.face_nodes
reference_fe_interior_nodes(e::ReferenceFEStencil) = e.interior_nodes
reference_fe_type(::ReferenceFEStencil{IType, RType, RefFE}) where {IType, RType, RefFE} = RefFE
reference_fe_vertex_nodes(e::ReferenceFEStencil) = e.vertex_nodes

function Base.show(io::IO, e::E) where {E <: ReferenceFEStencil}
  str = "ParentElement:\n"
  str = str * "  Coordinates: "
  for coord in eachcol(e.coordinates)
    str = str * "$coord\n               "
  end
  str = str * "\n"
  str = str * "  Vertex nodes: $(e.vertex_nodes)\n\n"
  str = str * "  Face nodes: "
  for face_nodes in eachcol(e.face_nodes)
    str = str * "$face_nodes\n              "
  end
  str = str * "\n"
  str = str * "  Interior nodes: $(e.interior_nodes)\n" 
  print(io, str)
end
