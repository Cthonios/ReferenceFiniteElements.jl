module ReferenceFiniteElementsSymbolicsExt

using ReferenceFiniteElements
using StaticArrays
using StructArrays
using Symbolics

function symbolic_shape_function_gradients(e::ReferenceFiniteElements.AbstractElementType, X, ξ, backend)
  Ns = symbolic_shape_function_values(e, X, ξ, backend)
  Ns = map(x -> Symbolics.derivative(x, ξ), Ns)
  return ReferenceFiniteElements.convert_to_matrix(e, backend, Ns...)
end

function symbolic_shape_function_hessians(e::ReferenceFiniteElements.AbstractElementType, X, ξ, backend)
  Ns = symbolic_shape_function_values(e, X, ξ, backend)
  Ns = map(x -> Symbolics.derivative(Symbolics.derivative(x, ξ), ξ), Ns)
  return ReferenceFiniteElements.convert_to_matrix(e, backend, Ns...)
end

function symbolic_shape_function_values(e::Edge{Lagrange, P, Q}, X, ξ, backend) where {P, Q}
  Ns = Vector{Num}(undef, num_nodes(e))
  for n in axes(Ns, 1)
    Ns[n] = 1
    for i in axes(Ns, 1)
      if i == n
        continue
      end
      Ns[n] = Ns[n] * (ξ - X[i][1]) / (X[n][1] - X[i][1])
    end
  end

  return ReferenceFiniteElements.convert_to_vector(e, backend, Ns...)
end

# function symbolic_shape_function_gradients(e::Edge{Lagrange, P, Q}, X, ξ, backend) where {P, Q}
#   Ns = symbolic_shape_function_values(e, X, ξ, backend)
#   return map(x -> Symbolics.derivative(x, ξ), Ns)
# end

# function symbolic_shape_function_hessians(e::Edge{Lagrange, P, Q}, X, ξ, backend) where {P, Q}
#   Ns = symbolic_shape_function_values(e, X, ξ, backend)
#   return map(x -> Symbolics.derivative(Symbolics.derivative(x, ξ), ξ), Ns)
# end

function ReferenceFiniteElements.CellInterpolants{Num}(e::ReferenceFiniteElements.AbstractElementType, Xs, backend)
  @variables ξ
  ξs, ws = ReferenceFiniteElements.quadrature_points_and_weights(e, backend)
  Ns = symbolic_shape_function_values(e, Xs, ξ, backend)
  Ns = fill(Ns, num_quadrature_points(e))
  ∇N_ξs = symbolic_shape_function_gradients(e, Xs, ξ, backend)
  ∇N_ξs = fill(∇N_ξs, num_quadrature_points(e))
  ∇∇N_ξs = symbolic_shape_function_hessians(e, Xs, ξ, backend)
  ∇∇N_ξs = fill(∇∇N_ξs, num_quadrature_points(e))
  vals = StructArray{ReferenceFiniteElements.Interpolants}((ξs, ws, Ns, ∇N_ξs, ∇∇N_ξs))
  return ReferenceFiniteElements.CellInterpolants(vals)
end

function ReferenceFiniteElements.SurfaceInterpolants{Num}(e::ReferenceFiniteElements.AbstractElementType, Xs, backend)
  @variables ξ
  ξs, ws = ReferenceFiniteElements.surface_quadrature_points_and_weights(e, backend)
  ξs = mapreduce(x -> x, hcat, ξs)
  ws = mapreduce(x -> x, hcat, ws)
  Ns = symbolic_shape_function_values(e, Xs, ξ, backend)
  Ns = fill(Ns, size(ξs))
  ∇N_ξs = symbolic_shape_function_gradients(e, Xs, ξ, backend)
  ∇N_ξs = fill(∇N_ξs, size(ξs))
  ∇∇N_ξs = symbolic_shape_function_hessians(e, Xs, ξ, backend)
  ∇∇N_ξs = fill(∇∇N_ξs, size(ξs))
  vals = StructArray{ReferenceFiniteElements.Interpolants}((ξs, ws, Ns, ∇N_ξs, ∇∇N_ξs))
  return ReferenceFiniteElements.SurfaceInterpolants(vals)
end

function ReferenceFiniteElements.ReferenceFE{Itype, Ftype, T, Num}(e) where {Itype, Ftype, T, Num}
  surf_e = ReferenceFiniteElements.surface_element(e)
  backend = ReferenceFiniteElements.ArrayBackend{T}()
  edge_nodes = ReferenceFiniteElements.element_edge_nodes(e, backend)
  face_nodes = ReferenceFiniteElements.element_face_nodes(e, backend)
  interior_nodes = ReferenceFiniteElements.element_interior_nodes(e, backend)
  Xs = ReferenceFiniteElements.nodal_coordinates(e, backend)
  interps = ReferenceFiniteElements.CellInterpolants{Num}(e, Xs, backend)
  surface_Xs = ReferenceFiniteElements.surface_nodal_coordinates(e, backend)
  surface_interps = ReferenceFiniteElements.SurfaceInterpolants{Num}(e, Xs, backend)
  return ReferenceFE{
    Itype, Ftype, typeof(e), typeof(surf_e), typeof(backend),
    typeof(edge_nodes), typeof(face_nodes), typeof(interior_nodes),
    typeof(Xs), typeof(interps),
    typeof(surface_Xs), typeof(surface_interps)
  }(
    e, surf_e, backend, 
    edge_nodes, face_nodes, interior_nodes, 
    Xs, interps,
    surface_Xs, surface_interps
  )
end

ReferenceFiniteElements.ReferenceFE{E, I, P, Q, Num}() where {E, I, P, Q} = 
ReferenceFE{Int64, Float64, SArray, Num}(E{I, P, Q}())

end # module