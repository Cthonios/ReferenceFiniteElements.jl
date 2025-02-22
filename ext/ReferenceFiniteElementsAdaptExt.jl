module ReferenceFiniteElementsAdaptExt

using Adapt
using ReferenceFiniteElements

function Adapt.adapt_structure(to, interp::ReferenceFiniteElements.Interpolants)
  ξ = Adapt.adapt_structure(to, interp.ξ)
  w = Adapt.adapt_structure(to, interp.w)
  N = Adapt.adapt_structure(to, interp.N)
  ∇N_ξ = Adapt.adapt_structure(to, interp.∇N_ξ)
  ∇∇N_ξ = Adapt.adapt_structure(to, interp.∇∇N_ξ)
  return ReferenceFiniteElements.Interpolants(ξ, w, N, ∇N_ξ, ∇∇N_ξ)
end

function Adapt.adapt_structure(to, interp::ReferenceFiniteElements.CellInterpolants)
  vals = Adapt.adapt_structure(to, interp.vals)
  return ReferenceFiniteElements.CellInterpolants(vals)
end

function Adapt.adapt_structure(to, interp::ReferenceFiniteElements.SurfaceInterpolants)
  vals = Adapt.adapt_structure(to, interp.vals)
  return ReferenceFiniteElements.SurfaceInterpolants(vals)
end

function Adapt.adapt_structure(to, re::ReferenceFE)
  I = ReferenceFiniteElements.integer_type(re)
  F = ReferenceFiniteElements.float_type(re)

  element = re.element
  surface_element = re.surface_element
  backend = re.backend

  edge_nodes = Adapt.adapt_structure(to, re.edge_nodes)
  face_nodes = Adapt.adapt_structure(to, re.face_nodes)
  interior_nodes = Adapt.adapt_structure(to, re.interior_nodes)
  Xs = Adapt.adapt_structure(to, re.Xs)
  cell_interps = Adapt.adapt_structure(to, re.cell_interps)
  surface_Xs = Adapt.adapt_structure(to, re.surface_Xs)
  surface_interps = Adapt.adapt_structure(to, re.surface_interps)

  return ReferenceFE{
    I, F, typeof(element), typeof(surface_element), typeof(backend),
    typeof(edge_nodes), typeof(face_nodes), typeof(interior_nodes),
    typeof(Xs), typeof(cell_interps),
    typeof(surface_Xs), typeof(surface_interps)
  }(
    element, surface_element, backend, 
    edge_nodes, face_nodes, interior_nodes, 
    Xs, cell_interps,
    surface_Xs, surface_interps
  )
end

end # module
