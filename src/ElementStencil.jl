struct ElementStencil{IType <: Integer, RType <: Real, RefFE <: AbstractReferenceFE}
  element_type::RefFE
  degree::IType
  coordinates::VecOrMat{RType}
  vertex_nodes::Vector{IType}
  face_nodes::Matrix{IType}
  interior_nodes::Vector{IType}
end
function Base.show(io::IO, e::E) where {E <: ElementStencil}
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

# @recipe function f(e::E) where E <: ElementStencil
#   # legend --> "off"

# end

export ElementStencil
