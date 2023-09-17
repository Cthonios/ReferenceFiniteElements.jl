module ReferenceFiniteElementsAdaptExt

using Adapt
using ReferenceFiniteElements

function Adapt.adapt_structure(
  to, 
  re::ReferenceFiniteElements.ReferenceFE{Itype, N, D, Ftype, L1, L2, S}
) where {Itype, N, D, Ftype, L1, L2, S}

  nodal_coordinates = Adapt.adapt_structure(to, re.nodal_coordinates)
  face_nodes        = Adapt.adapt_structure(to, re.face_nodes)
  interior_nodes    = Adapt.adapt_structure(to, re.interior_nodes)
  interpolants      = Adapt.adapt_structure(to, re.interpolants)

  ReferenceFiniteElements.ReferenceFE{Itype, N, D, Ftype, L1, L2, typeof(interpolants),
                                      typeof(nodal_coordinates), typeof(face_nodes), typeof(interior_nodes)}(
    nodal_coordinates, face_nodes, interior_nodes, interpolants
  )
end

end # module