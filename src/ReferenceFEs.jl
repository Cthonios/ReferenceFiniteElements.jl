# main type of the package
"""
ReferenceFE is the main type. 

This type defines the nodal coordinates, face nodes, interior nodes
and interpolants
"""
# struct ReferenceFE{
#   Itype, N, D, Ftype, L1, L2, S, 
#   RefFE <: ReferenceFEType{N, D},
#   VOM <: AbstractVecOrMat, 
#   M <: AbstractMatrix, 
#   V <: AbstractVector
# }
struct ReferenceFE{
  Itype, N, D, Ftype, L1, L2, Q,
  RefFEType <: ReferenceFEType{N, D},
  S,
  VOM       <: AbstractVecOrMat,
  M         <: AbstractMatrix,
  V         <: AbstractVector
}
  ref_fe_type::RefFEType
  nodal_coordinates::VOM
  face_nodes::M
  interior_nodes::V
  interpolants::S
end

"""
Constructor for ReferenceFE
"""
function ReferenceFE(
  e::ReferenceFEType{N, D, Q};
  int_type::Type{<:Integer} = Int64, 
  float_type::Type{<:Number} = Float64,
  array_type::Type{<:Union{<:MArray, <:SArray}} = SArray
) where {N, D, Q}

  nodal_coordinates, face_nodes, interior_nodes = element_stencil(e, int_type, float_type)
  interps = Interpolants(e, array_type, float_type)

  return ReferenceFE{
    int_type, N, D, float_type, N * D, N * D * D, Q,
    typeof(e), typeof(interps), typeof(nodal_coordinates),
    typeof(face_nodes), typeof(interior_nodes)
  }(e, nodal_coordinates, face_nodes, interior_nodes, interps)
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
element_type(::ReferenceFE{
  Itype, N, D, Ftype, L1, L2, Q, RefFE, S, VOM, M, V
}) where {Itype, N, D, Ftype, L1, L2, Q, RefFE <: ReferenceFEType, S, VOM, M, V} = RefFE

vertex_nodes(::ReferenceFE{
  Itype, N, D, Ftype, L1, L2, Q, RefFE, S, VOM, M, V
}) where {Itype, N, D, Ftype, Q, L1, L2, RefFE <: ReferenceFEType, S, VOM, M, V} = Base.OneTo(N)

"""
Returns number of nodes per element
"""
num_nodes_per_element(e::ReferenceFE) = num_nodes(e.ref_fe_type)

"""
Returns number of quadrature points
"""
num_q_points(e::ReferenceFE) = num_q_points(e.ref_fe_type)