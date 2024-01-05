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
# struct ReferenceFE{
#   Itype, N, D, Ftype, L1, L2, Q,
#   RefFEType <: ReferenceFEType{N, D},
#   S,
#   VOM       <: AbstractVecOrMat,
#   M         <: AbstractMatrix,
#   V         <: AbstractVector
# }
struct ReferenceFE{
  Itype, Ftype, N, D, Q,
  RefFEType <: ReferenceFEType{N, D},
  S         <: StructArray, # TODO figure out how to genarilze beyond structarrays
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
  array_type::Type{<:Union{<:MArray, <:SArray, <:Array}} = SArray
) where {N, D, Q}

  nodal_coordinates, face_nodes, interior_nodes = element_stencil(e, int_type, float_type)
  interps = Interpolants{array_type, float_type}(e)

  return ReferenceFE{
    int_type, float_type, N, D, Q,
    typeof(e), typeof(interps), typeof(nodal_coordinates),
    typeof(face_nodes), typeof(interior_nodes)
  }(e, nodal_coordinates, face_nodes, interior_nodes, interps)
end

"""
"""
function Base.show(io::IO, e::ReferenceFE)
  print(io, "ReferenceFE\n")
  print(io, "  Element type                = $(e.ref_fe_type |> typeof)\n")
  print(io, "  Dimension                   = $(num_dimensions(e))\n")
  print(io, "  Number of nodes             = $(num_nodes_per_element(e))\n")
  print(io, "  Number of quadrature points = $(num_q_points(e))\n")
  print(io, "  Integer type                = $(int_type(e))\n")
  print(io, "  Float type                  = $(float_type(e))\n")
  print(io, "  Nodal coordinates type      = $(typeof(e.nodal_coordinates))\n")
  print(io, "  Face nodes type             = $(typeof(e.face_nodes))\n")
  print(io, "  Interior nodes type         = $(typeof(e.interior_nodes))\n")
  print(io, "  Interpolants type           = $(typeof(e.interpolants).name.name){$(typeof(e.interpolants[1]).name.name)}\n")
  print(io, "\n")
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
# element_type(::ReferenceFE{
#   Itype, N, D, Ftype, L1, L2, Q, RefFE, S, VOM, M, V
# }) where {Itype, N, D, Ftype, L1, L2, Q, RefFE <: ReferenceFEType, S, VOM, M, V} = RefFE

# vertex_nodes(::ReferenceFE{
#   Itype, N, D, Ftype, L1, L2, Q, RefFE, S, VOM, M, V
# }) where {Itype, N, D, Ftype, Q, L1, L2, RefFE <: ReferenceFEType, S, VOM, M, V} = Base.OneTo(N)

# int_type(::ReferenceFE{
#   Itype, N, D, Ftype, L1, L2, Q, RefFE, S, VOM, M, V
# }) where {Itype, N, D, Ftype, Q, L1, L2, RefFE <: ReferenceFEType, S, VOM, M, V} = Itype

# float_type(::ReferenceFE{
#   Itype, N, D, Ftype, L1, L2, Q, RefFE, S, VOM, M, V
# }) where {Itype, N, D, Ftype, Q, L1, L2, RefFE <: ReferenceFEType, S, VOM, M, V} = Ftype

element_type(::ReferenceFE{
  Itype, Ftype, N, D, Q, RefFE, S, VOM, M, V
}) where {Itype, Ftype, N, D, Q, RefFE <: ReferenceFEType, S, VOM, M, V} = RefFE

vertex_nodes(::ReferenceFE{
  Itype, Ftype, N, D, Q, RefFE, S, VOM, M, V
}) where {Itype, Ftype, N, D, Q, RefFE <: ReferenceFEType, S, VOM, M, V} = Base.OneTo(N)

int_type(::ReferenceFE{
  Itype, Ftype, N, D, Q, RefFE, S, VOM, M, V
}) where {Itype, Ftype, N, D, Q, RefFE <: ReferenceFEType, S, VOM, M, V} = Itype

float_type(::ReferenceFE{
  Itype, Ftype, N, D, Q, RefFE, S, VOM, M, V
}) where {Itype, Ftype, N, D, Q, RefFE <: ReferenceFEType, S, VOM, M, V} = Ftype


"""
Returns number of dimensions
"""
num_dimensions(e::ReferenceFE) = num_dimensions(e.ref_fe_type)

"""
Returns number of nodes per element
"""
num_nodes_per_element(e::ReferenceFE) = num_nodes(e.ref_fe_type)

"""
Returns number of quadrature points
"""
num_q_points(e::ReferenceFE) = num_q_points(e.ref_fe_type)

