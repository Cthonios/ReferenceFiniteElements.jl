module ReferenceFiniteElementsMakieCoreExt

using MakieCore
using ReferenceFiniteElements

"""
This method is for plotting the vertices
"""
function MakieCore.convert_arguments(
  p::Type{<:Scatter},
  re::ReferenceFE{Itype, Ftype, N, 2, Q, RefFEType, Interp, Vom, M, V}
) where {Itype, Ftype, N, Q, RefFEType, Interp, Vom, M, V}
  nodal_coords = re.nodal_coordinates
  return MakieCore.convert_arguments(p, nodal_coords[1, :], nodal_coords[2, :])
end 

function MakieCore.convert_arguments(
  p::Type{<:LineSegments},
  re::ReferenceFE{Itype, Ftype, N, 2, Q, RefFEType, Interp, Vom, M, V}
) where {Itype, Ftype, N, Q, RefFEType, Interp, Vom, M, V}
  nodal_coords = re.nodal_coordinates
  edge_nodes = re.edge_nodes

  return MakieCore.convert_arguments(p, nodal_coords[1, :], nodal_coords[2, :])
end

end # module