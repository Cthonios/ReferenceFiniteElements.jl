struct ReferenceFE{
  Itype, Ftype, Etype, SEType, Backend,
  EdgeNodes, FaceNodes, InteriorNodes,
  Coords, CellInterps,
  SurfaceCoords, SurfaceInterps
}
  element::Etype
  surface_element::SEType
  backend::Backend
  edge_nodes::EdgeNodes
  face_nodes::FaceNodes
  interior_nodes::InteriorNodes
  Xs::Coords
  cell_interps::CellInterps
  surface_Xs::SurfaceCoords
  surface_interps::SurfaceInterps
end

function ReferenceFE{Itype, Ftype, T}(e::AbstractElementType) where {Itype, Ftype, T}
  surf_e = surface_element(e)
  backend = ArrayBackend{T}()
  edge_nodes = element_edge_nodes(e, backend)
  face_nodes = element_face_nodes(e, backend)
  interior_nodes = element_interior_nodes(e, backend)
  Xs = nodal_coordinates(e, backend)
  interps = CellInterpolants(e, Xs, backend)
  surface_Xs = surface_nodal_coordinates(e, backend)
  surface_interps = SurfaceInterpolants(e, Xs, backend)

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

ReferenceFE(e::AbstractElementType) = ReferenceFE{Int64, Float64, SArray}(e)
ReferenceFE{E, I, P, Q}() where {E, I, P, Q} = ReferenceFE{Int64, Float64, SArray}(E{I, P, Q}())
ReferenceFE{T, E, I, P, Q}() where {T, E, I, P, Q} = ReferenceFE{Int64, Float64, T}(E{I, P, Q}())
ReferenceFE{T1, T2, T3, E, I, P, Q}() where {T1, T2, T3, E, I, P, Q} = ReferenceFE{T1, T2, T3}(E{I, P, Q}())

function Base.show(io::IO, re::ReferenceFE; tab="")
  println(io, "$(tab)ReferenceFE:")
  println(io, "$(tab)  Element Type                                  = $(typeof(re.element))")
  println(io, "$(tab)  Surface Element Type                          = $(typeof(re.surface_element))")
  println(io, "$(tab)  Interpolation Type                            = $(interpolation_type(re.element))") 
  println(io, "$(tab)  Polynomial Degree                             = $(polynomial_degree(re.element))")
  println(io, "$(tab)  Quadrature Degree                             = $(quadrature_degree(re.element))")
  println(io, "$(tab)  Array Backend Type                            = $(re.backend)")
  println(io, "$(tab)  Edge Nodes Storage Type                       = $(typeof(re.edge_nodes))")
  println(io, "$(tab)  Face Nodes Storage Type                       = $(typeof(re.face_nodes))")
  println(io, "$(tab)  Interior Nodes Storage Type                   = $(typeof(re.interior_nodes))")
  println(io, "$(tab)  Nodal Coordinates Storage Type                = $(typeof(re.Xs))")
  show(io, re.cell_interps; tab="  ")
  println(io, "$(tab)  Surface Nodal Coordinates Storage Type        = $(typeof(re.surface_Xs))")
  show(io, re.surface_interps; tab="  ")
end

integer_type(::ReferenceFE{I, F}) where {I, F} = I
float_type(::ReferenceFE{I, F}) where {I, F} = F

# size helpers
dimension(e::ReferenceFE) = dimension(e.element)
element_type(::ReferenceFE{Itype, Ftype, Etype}) where {Itype, Ftype, Etype} = Etype
interpolation_type(e::ReferenceFE) = interpolation_type(e::ReferenceFE)
num_nodes(e::ReferenceFE) = num_nodes(e.element)
num_shape_functions(e::ReferenceFE) = num_shape_functions(e.element)
num_quadrature_points(e::ReferenceFE) = length(e.cell_interps.vals.w)
polynomial_degree(e::ReferenceFE) = polynomial_degree(e.element)
quadrature_degree(e::ReferenceFE) = quadrature_degree(e.element)

# getters
quadrature_point(e::ReferenceFE, q) = e.cell_interps.vals.ξ[q]
quadrature_points(e::ReferenceFE) = e.cell_interps.vals.ξ

quadrature_weight(e::ReferenceFE, q) = e.cell_interps.vals.w[q]
quadrature_weights(e::ReferenceFE) = e.cell_interps.vals.w

shape_function_gradient(e::ReferenceFE, q) = e.cell_interps.vals.∇N_ξ[q]
shape_function_gradients(e::ReferenceFE) = e.cell_interps.vals.∇N_ξ

shape_function_hessian(e::ReferenceFE, q) = e.cell_interps.vals.∇∇N_ξ[q]
shape_function_hessians(e::ReferenceFE) = e.cell_interps.vals.∇∇N_ξ

shape_function_value(e::ReferenceFE, q) = e.cell_interps.vals.N[q]
shape_function_values(e::ReferenceFE) = e.cell_interps.vals.N

surface_quadrature_point(e::ReferenceFE, q, f) = e.surface_interps.vals.ξ[q, f]
surface_quadrature_points(e::ReferenceFE) = e.surface_interps.vals.ξ

surface_quadrature_weight(e::ReferenceFE, q, f) = e.surface_interps.vals.w[q, f]
surface_quadrature_weights(e::ReferenceFE) = e.surface_interps.vals.w

surface_shape_function_gradient(e::ReferenceFE, q, f) = e.surface_interps.vals.∇N_ξ[q, f]
surface_shape_function_gradients(e::ReferenceFE) = e.surface_interps.vals.∇N_ξ

surface_shape_function_hessian(e::ReferenceFE, q, f) = e.surface_interps.vals.∇∇N_ξ[q, f]
surface_shape_function_hessians(e::ReferenceFE) = e.surface_interps.vals.∇∇N_ξ

surface_shape_function_value(e::ReferenceFE, q, f) = e.surface_interps.vals.N[q, f]
surface_shape_function_values(e::ReferenceFE) = e.surface_interps.vals.N

surface_normal(e::ReferenceFE, q, f) = surface_quadrature_point(e, q, f)
surface_normals(e::ReferenceFE, q, f) = surface_quadrature_point(e, q, f)

# other methods for working in physical space
function mapping_jacobian(e::ReferenceFE, X, q)
  J = (X' * shape_function_gradient(e, q))'
  return J
end

function calculate_JxW(e::ReferenceFE, X, q)
  J = mapping_jacobian(e, X, q)
  w = quadrature_weight(e, q)
  return w * det(J)
end

function map_shape_function_gradient(e::ReferenceFE, X, q)
  ∇N_ξ = shape_function_gradient(e, q)
  J = (X' * ∇N_ξ)'
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'
  return ∇N_X
end

# function map_shape_function_gradient(e::ReferenceFE, interp, X)
#   return Interpolants(e.ξ, e.w, e.N, map_shape_function_gradient(e, X, interp.∇N_ξ))
# end

# below assumes integrand values are already calculated
function integrate(e::ReferenceFE, vals)
  return sum(quadrature_weights(e) .* vals)
end

# TODO need more general integration method
# that can handle a function and take in maybe X
# as an input with general number of args supported
# maybe something like
# f(X, N, ∇N, ∇∇N, args...)

function integrate(f, e::ReferenceFE, X, args...)
  val = zero(eltype(X))
  for q in 1:num_quadrature_points(e)
    val = val + calculate_JxW(e, X, q) * f(
      X, 
      shape_function_value(e, q), 
      map_shape_function_gradient(e, X, q), 
      args...
    )
  end
  return val
end
