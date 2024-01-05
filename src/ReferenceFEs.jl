# main type of the package
"""
ReferenceFE is the main type. 

This type defines the nodal coordinates, face nodes, interior nodes
and interpolants
"""
struct ReferenceFE{
  Itype, Ftype, N, D, Q,
  RefFEType <: ReferenceFEType{N, D},
  Interp,
  VOM       <: AbstractVecOrMat{Ftype},
  # M1        <: AbstractArray, # can be an empty array
  # M2        <: AbstractArray,
  M1        <: AbstractArray{Itype, 2},
  M2        <: AbstractArray{Itype, 2},
  V         <: AbstractArray{Itype, 1}
}
  ref_fe_type::RefFEType
  nodal_coordinates::VOM
  edge_nodes::M1
  face_nodes::M2
  interior_nodes::V
  interpolants::Interp
end

"""
Constructor for ReferenceFE
"""
function ReferenceFE(
  e::ReferenceFEType{N, D, Q};
  int_type::Type{<:Integer} = Int64, 
  float_type::Type{<:Number} = Float64,
  array_type::Type{<:Union{<:MArray, <:SArray, <:Array}} = SArray,
  storage_type::Type{<:Union{<:Array, <:StructArray}} = StructArray
) where {N, D, Q}

  nodal_coordinates, edge_nodes, face_nodes, interior_nodes = element_stencil(e, int_type, float_type)
  interps = Interpolants{storage_type, array_type, float_type}(e)

  return ReferenceFE{
    int_type, float_type, N, D, Q,
    typeof(e), typeof(interps), typeof(nodal_coordinates),
    typeof(edge_nodes), typeof(face_nodes), typeof(interior_nodes)
  }(e, nodal_coordinates, edge_nodes, face_nodes, interior_nodes, interps)
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
  print(io, "  Edge nodes type             = $(typeof(e.edge_nodes))\n")
  print(io, "  Face nodes type             = $(typeof(e.face_nodes))\n")
  print(io, "  Interior nodes type         = $(typeof(e.interior_nodes))\n")
if isa(e.interpolants, Interpolants)
    print(io, "  Interpolants type           = Interpolants{$(eltype(e.interpolants.N).name.name), ",
              "$(eltype(e.interpolants.∇N_ξ).name.name), ",
              "$(eltype(e.interpolants.∇∇N_ξ).name.name)}"
    )
  else
    print(io, "  Interpolants type           = $(typeof(e.interpolants).name.name){Interpolants}\n")
  end
  print(io, "\n")
end

"""
Returns all quadrature points
"""
quadrature_points(e::ReferenceFE) = e.interpolants.ξ

"""
Returns a specific quadrature point indexed by q
"""
# quadrature_point(e::ReferenceFE, q::Int) = LazyRow(e.interpolants, q).ξ
quadrature_points(e::ReferenceFE, q::Int) = e.interpolants.ξ[q]

"""
Returns all quadrature weights
"""
quadrature_weights(e::ReferenceFE) = e.interpolants.w

"""
Returns a specific quadrature weight indexed by q
"""
# quadrature_weight(e::ReferenceFE, q::Int) = LazyRow(e.interpolants, q).w
quadrature_weights(e::ReferenceFE, q::Int) = e.interpolants.w[q]

"""
Returns all shape function values
"""
shape_function_values(e::ReferenceFE) = e.interpolants.N

"""
Returns a specific quadrature point's shape function value indexed by q
"""
shape_function_values(e::ReferenceFE, q::Int) = e.interpolants.N[q]

"""
Returns all shape function gradients
"""
shape_function_gradients(e::ReferenceFE) = e.interpolants.∇N_ξ

"""
Returns a specific quadrature point's shape function gradient indexed by q
"""
shape_function_gradients(e::ReferenceFE, q::Int) = e.interpolants.∇N_ξ[q]

"""
Returns all shape function hessians
"""
shape_function_hessians(e::ReferenceFE) = e.interpolants.∇∇N_ξ

"""
Returns a specific quadrature point's shape function hessian indexed by q
"""
shape_function_hessians(e::ReferenceFE, i::Int) = e.interpolants.∇∇N_ξ[i]


"""
Returns the element type
"""
element_type(::ReferenceFE{
  Itype, Ftype, N, D, Q, RefFE, S, VOM, M1, M2, V
}) where {Itype, Ftype, N, D, Q, RefFE <: ReferenceFEType, S, VOM, M1, M2, V} = RefFE

"""
Returns the nodes of vertices
TODO this is probably not very useful
"""
vertex_nodes(::ReferenceFE{
  Itype, Ftype, N, D, Q, RefFE, S, VOM, M1, M2, V
}) where {Itype, Ftype, N, D, Q, RefFE <: ReferenceFEType, S, VOM, M1, M2, V} = Base.OneTo(N)

"""
Returns the integer type used to store node ids and such
"""
int_type(::ReferenceFE{
  Itype, Ftype, N, D, Q, RefFE, S, VOM, M1, M2, V
}) where {Itype, Ftype, N, D, Q, RefFE <: ReferenceFEType, S, VOM, M1, M2, V} = Itype

"""
Returns the float type used to store nodal coordinates and interpolation arrays
"""
float_type(::ReferenceFE{
  Itype, Ftype, N, D, Q, RefFE, S, VOM, M1, M2, V
}) where {Itype, Ftype, N, D, Q, RefFE <: ReferenceFEType, S, VOM, M1, M2, V} = Ftype


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

"""
Returns the edge list
"""
edges(e::ReferenceFE) = e.edges