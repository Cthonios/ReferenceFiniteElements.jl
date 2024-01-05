module ReferenceFiniteElementsAdaptExt

using Adapt
using ReferenceFiniteElements

function Adapt.adapt_structure(to, re::ReferenceFE)
  Itype = ReferenceFiniteElements.int_type(re)
  N     = num_nodes_per_element(re)
  D     = num_dimensions(re)
  Ftype = ReferenceFiniteElements.float_type(re)
  Q     = num_q_points(re)
  
  ref_fe_type       = Adapt.adapt_structure(to, re.ref_fe_type)
  nodal_coordinates = Adapt.adapt_structure(to, re.nodal_coordinates)
  edge_nodes        = Adapt.adapt_structure(to, re.edge_nodes)
  face_nodes        = Adapt.adapt_structure(to, re.face_nodes)
  interior_nodes    = Adapt.adapt_structure(to, re.interior_nodes)
  interpolants      = Adapt.adapt_structure(to, re.interpolants)

  return ReferenceFE{
    # Itype, N, D, Ftype, N * D, N * D * D, Q, 
    Itype, Ftype, N, D, Q,
    typeof(ref_fe_type), typeof(interpolants), typeof(nodal_coordinates), 
    typeof(edge_nodes), typeof(face_nodes), typeof(interior_nodes)
  }(ref_fe_type, nodal_coordinates, edge_nodes, face_nodes, interior_nodes, interpolants)
end

end # module