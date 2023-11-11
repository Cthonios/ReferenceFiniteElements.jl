module ReferenceFiniteElements

# element types
# export Edge
export Hex8,
       Quad4, Quad9,
       Tet4, Tet10,
       Tri3, Tri6

# types
export ReferenceFE

# methods
export quadrature_point,
       quadrature_points,
       quadrature_weight,
       quadrature_weights,
       shape_function_gradients,
       shape_function_hessians,
       shape_function_values,
       vertex_nodes

# dependencies
using DocStringExtensions,
      FastGaussQuadrature,
      LinearAlgebra,
      Polynomials,
      PrecompileTools,
      SpecialPolynomials,
      StaticArrays,
      StructArrays

# for docs
@template (FUNCTIONS, METHODS, MACROS) =
"""
$(TYPEDSIGNATURES)
$(DOCSTRING)
$(METHODLIST)
"""

@template (TYPES) =
"""
$(TYPEDFIELDS)
$(DOCSTRING)
"""

# type used to exploit multiple dispatch
"""
Type to define new element shapes
"""
abstract type ReferenceFEType{N, D} end

"""
Returns the quadrature degree of a ReferenceFEType
"""
degree(e::ReferenceFEType) = e.degree

"""
Returns the number of nodes for a ReferenceFEType
"""
num_nodes(::ReferenceFEType{N, D}) where {N, D} = N

"""
Returns the number of dimensions for a ReferenceFEType
"""
num_dimensions(::ReferenceFEType{N, D}) where {N, D} = D

# type to aid in making ReferenceFE with StructArrays.jl
"""
Interpolant container for a single quadrature point
"""
struct Interpolants{N, D, Ftype, L1, L2}
  ξ::SVector{D, Ftype}
  w::Ftype
  N::SVector{N, Ftype}
  ∇N_ξ::SMatrix{N, D, Ftype, L1}
  ∇∇N_ξ::SArray{Tuple{N, D, D}, Ftype, 3, L2}
end

"""
Constructor for a StructArray of Interpolants
TODO maybe change the name of this. It might be confusing
"""
function Interpolants(
  e::ReferenceFEType{N, D}, ::Type{Ftype} = Float64
) where {N, D, Ftype}

  ξs_temp, ws = quadrature_points_and_weights(e, Ftype)
  ξs = reinterpret(SVector{D, Ftype}, vec(ξs_temp))
  Ns = Vector{SVector{N, Ftype}}(undef, length(ξs))
  ∇N_ξs = Vector{SMatrix{N, D, Ftype, N * D}}(undef, length(ξs))
  ∇∇N_ξs = Vector{SArray{Tuple{N, D, D}, Ftype, 3, N * D * D}}(undef, length(ξs))
  for (q, ξ) in enumerate(ξs)
    Ns[q]     = shape_function_values(e, ξ)
    ∇N_ξs[q]  = shape_function_gradients(e, ξ)
    ∇∇N_ξs[q] = shape_function_hessians(e, ξ)
  end
  s = StructArray{Interpolants{N, D, Ftype, N * D, N * D * D}}((
    ξ=ξs, w=ws, N=Ns, ∇N_ξ=∇N_ξs, ∇∇N_ξ=∇∇N_ξs
  ))
  return s
end

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
  print(io, "Nodal coordinates        = \n")
  display(e.nodal_coordinates)
  print(io, "Face nodes               = \n")
  display(e.face_nodes)
  print(io, "\n")
  display(e.interior_nodes)
  print(io, "\n")
  print(io, "Shape function values    = \n")
  Ns = shape_function_values(e)
  for n in axes(e.interpolants, 1)
    # display(shape_function_values(e, n))
    display(Ns[n])
  end
  print(io, "\n")
  print(io, "Shape function gradients = \n")
  ∇N_ξs = shape_function_gradients(e)
  for n in axes(e.interpolants, 1)
    # display(shape_function_gradients(e, n))
    display(∇N_ξs[n])
  end
  print(io, "\n")
  print(io, "Shape function hessians  = \n")
  ∇∇N_ξs = shape_function_hessians(e)
  for n in axes(e.interpolants, 1)
    # display(shape_function_hessians(e, n))
    display(∇∇N_ξs[n])
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

# implementations of things common across multiple element types
include("implementations/HexCommon.jl")
include("implementations/QuadCommon.jl")
include("implementations/TetCommon.jl")
include("implementations/TriCommon.jl")

# implementations of things specific to element types
# include("implementations/Edge.jl")
include("implementations/Hex8.jl")
include("implementations/Quad4.jl")
include("implementations/Quad9.jl")
include("implementations/Tet4.jl")
include("implementations/Tet10.jl")
include("implementations/Tri3.jl")
include("implementations/Tri6.jl")

# include("implementations/SimplexTri.jl")

# precompilation
@setup_workload begin
  @compile_workload begin
    for degree in [1, 2, 3, 4, 5, 6]
      ReferenceFE(Hex8(degree))
      ReferenceFE(Quad4(degree))
      ReferenceFE(Quad9(degree))
      ReferenceFE(Tri3(degree))
      ReferenceFE(Tri6(degree))
    end

    for degree in [1, 2]
      ReferenceFE(Tet4(degree))
      ReferenceFE(Tet10(degree))
    end
  end
end

end # module ReferenceFEs
