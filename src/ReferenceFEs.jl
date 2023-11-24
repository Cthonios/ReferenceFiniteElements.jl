# main type of the package
"""
ReferenceFE is the main type. 

This type defines the nodal coordinates, face nodes, interior nodes
and interpolants
"""
struct ReferenceFE{Itype, N, D, Ftype, L1, L2, S, VOM <: AbstractVecOrMat, M <: AbstractMatrix, V <: AbstractVector}
  nodal_coordinates::VOM
  face_nodes::M
  interior_nodes::V
  interpolants::S
end

"""
Constructor for ReferenceFE
"""
function ReferenceFE(
  e::ReferenceFEType{N, D},
  ::Type{Itype} = Int64, ::Type{Ftype} = Float64
) where {N, D, Itype, Ftype}

  nodal_coordinates, face_nodes, interior_nodes = element_stencil(e, Itype, Ftype)
  interps = Interpolants(e, Ftype)

  return ReferenceFE{Itype, N, D, Ftype, N * D, N * D *D, typeof(interps), 
                     typeof(nodal_coordinates), typeof(face_nodes), typeof(interior_nodes)}(
    nodal_coordinates, face_nodes, interior_nodes,
    interps
  )
end

"""
"""
function Base.show(io::IO, e::ReferenceFE)
  print(io, "Element type             = $(typeof(e))\n\n")
  # print(io, "Nodal coordinates        = \n")
  # display(e.nodal_coordinates)
  # print(io, "Face nodes               = \n")
  # display(e.face_nodes)
  # print(io, "\n")
  # display(e.interior_nodes)
  # print(io, "\n")
  # print(io, "Shape function values    = \n")
  # Ns = shape_function_values(e)
  # for n in axes(e.interpolants, 1)
  #   # display(shape_function_values(e, n))
  #   display(Ns[n])
  # end
  # print(io, "\n")
  # print(io, "Shape function gradients = \n")
  # ∇N_ξs = shape_function_gradients(e)
  # for n in axes(e.interpolants, 1)
  #   # display(shape_function_gradients(e, n))
  #   display(∇N_ξs[n])
  # end
  # print(io, "\n")
  # print(io, "Shape function hessians  = \n")
  # ∇∇N_ξs = shape_function_hessians(e)
  # for n in axes(e.interpolants, 1)
  #   # display(shape_function_hessians(e, n))
  #   display(∇∇N_ξs[n])
  # end
  # print(io, "\n")
end

"""
Returns all quadrature points
"""
quadrature_points(e::ReferenceFE) = e.interpolants.ξ

"""
Returns a specific quadrature point indexed by q
"""
quadrature_point(e::ReferenceFE, q::Int) = LazyRow(e.interpolants, q).ξ

"""
Returns all quadrature weights
"""
quadrature_weights(e::ReferenceFE) = e.interpolants.w

"""
Returns a specific quadrature weight indexed by q
"""
quadrature_weight(e::ReferenceFE, q::Int) = LazyRow(e.interpolants, q).w

"""
Returns all shape function values
"""
shape_function_values(e::ReferenceFE) = e.interpolants.N

"""
Returns a specific quadrature point's shape function value indexed by q
"""
shape_function_values(e::ReferenceFE, i::Int) = LazyRow(e.interpolants, i).N

"""
Returns all shape function gradients
"""
shape_function_gradients(e::ReferenceFE) = e.interpolants.∇N_ξ

"""
Returns a specific quadrature point's shape function gradient indexed by q
"""
shape_function_gradients(e::ReferenceFE, i::Int) = LazyRow(e.interpolants, i).∇N_ξ

"""
Returns all shape function hessians
"""
shape_function_hessians(e::ReferenceFE) = e.interpolants.∇∇N_ξ

"""
Returns a specific quadrature point's shape function hessian indexed by q
"""
shape_function_hessians(e::ReferenceFE, i::Integer) = LazyRow(getfield(e, :interpolants), i).∇∇N_ξ

"""
Returns the nodes of vertices
TODO this is probably not very useful
"""
vertex_nodes(::ReferenceFE{Itype, N, D, Ftype, L1, L2}) where {Itype, N, D, Ftype, L1, L2} = 1:N
