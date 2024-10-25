module ReferenceFiniteElementsAdaptExt

using Adapt
using ReferenceFiniteElements
using StructArrays

function Adapt.adapt_structure(to, interp::ReferenceFiniteElements.Interpolants)
  ξ = Adapt.adapt_structure(to, interp.ξ)
  w = Adapt.adapt_structure(to, interp.w)
  N = Adapt.adapt_structure(to, interp.N)
  ∇N_ξ = Adapt.adapt_structure(to, interp.∇N_ξ)
  ∇∇N_ξ = Adapt.adapt_structure(to, interp.∇∇N_ξ)
  return ReferenceFiniteElements.Interpolants(ξ, w, N, ∇N_ξ, ∇∇N_ξ)
end

function Adapt.adapt_structure(to, interp::ReferenceFiniteElements.CellInterpolants)
  vals = replace_storage(to, interp.vals)
  return return ReferenceFiniteElements.CellInterpolants(vals)
end

function Adapt.adapt_structure(to, interp::ReferenceFiniteElements.SurfaceInterpolants)
  vals = replace_storage(to, interp.vals)
  return return ReferenceFiniteElements.SurfaceInterpolants(vals)
end

function Adapt.adapt_structure(to, re::ReferenceFE)
  I = ReferenceFiniteElements.integer_type(re)
  F = ReferenceFiniteElements.float_type(re)

  element = re.element
  surface_element = re.surface_element
  backend = ReferenceFiniteElements.ArrayBackend{to}()

  edge_nodes = Adapt.adapt_structure(to, re.edge_nodes)
  face_nodes = Adapt.adapt_structure(to, re.face_nodes)
  interior_nodes = Adapt.adapt_structure(to, re.interior_nodes)
  Xs = Adapt.adapt_structure(to, re.Xs)
  cell_interps = Adapt.adapt_structure(to, re.cell_interps)
  # ξs = Adapt.adapt_structure(to, re.ξs)
  # ws = Adapt.adapt_structure(to, re.ws)
  # Ns = Adapt.adapt_structure(to, re.Ns)
  # ∇N_ξs = Adapt.adapt_structure(to, re.∇N_ξs)
  # ∇∇N_ξs = Adapt.adapt_structure(to, re.∇∇N_ξs)
  # surface_Xs = Adapt.adapt_structure(to, re.surface_Xs)
  surface_Xs = Adapt.adapt_structure(to, re.surface_Xs)
  surface_interps = Adapt.adapt_structure(to, re.surface_interps)
  # surface_ξs = Adapt.adapt_structure(to, re.surface_ξs)
  # surface_ξs = re.surface_ξs
  # # surface_ws = Adapt.adapt_structure(to, re.surface_ws)
  # surface_ws = re.surface_ws
  # # surface_Ns = Adapt.adapt_structure(to, re.surface_Ns)
  # surface_Ns = re.surface_Ns
  # # surface_∇N_ξs = Adapt.adapt_structure(to, re.surface_∇N_ξs)
  # surface_∇N_ξs = re.surface_∇N_ξs
  # # surface_∇∇N_ξs = Adapt.adapt_structure(to, re.surface_∇∇N_ξs)
  # surface_∇∇N_ξs = re.surface_∇∇N_ξs

  return ReferenceFE{
    I, F, typeof(element), typeof(surface_element), typeof(backend),
    typeof(edge_nodes), typeof(face_nodes), typeof(interior_nodes),
    typeof(Xs), typeof(cell_interps),
    # typeof(ξs), typeof(ws),
    # typeof(Ns), typeof(∇N_ξs), typeof(∇∇N_ξs),
    typeof(surface_Xs), typeof(surface_interps)
    # typeof(surface_ξs), typeof(surface_ws),
    # typeof(surface_Ns), typeof(surface_∇N_ξs), typeof(surface_∇∇N_ξs)
  }(
    element, surface_element, backend, 
    edge_nodes, face_nodes, interior_nodes, 
    Xs, cell_interps,
    # ξs, ws, 
    # Ns, ∇N_ξs, ∇∇N_ξs,
    surface_Xs, surface_interps
    # surface_ξs, surface_ws,
    # surface_Ns, surface_∇N_ξs, surface_∇∇N_ξs
  )
end

end # module