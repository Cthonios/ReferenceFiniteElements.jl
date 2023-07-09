module ReferenceFiniteElements

# element types
# export Edge
export Hex8
export Quad4, Quad9
# export Tet4, Tet10
export Tet4
export Tri3, Tri6

# types
export ReferenceFE

# methods
export quadrature_point
export quadrature_points
export quadrature_weight
export quadrature_weights
export shape_function_gradients
export shape_function_values
export vertex_nodes

# dependencies
using FastGaussQuadrature
using InteractiveUtils
using LinearAlgebra
using Polynomials
using PrecompileTools
using SpecialPolynomials
using StaticArrays

struct ReferenceFEType{N, D}
end
num_nodes(::ReferenceFEType{N, D}) where {N, D} = N
num_dimensions(::ReferenceFEType{N, D}) where {N, D} = D

# new type below
struct ReferenceFE{N, D, L, Itype, Ftype <: AbstractFloat}
  # element nodal info
  nodal_coordinates::VecOrMat{Ftype}
  # nodal_coordinates::Matrix{Ftype}
  face_nodes::Matrix{Itype}
  interior_nodes::Vector{Itype}
  # quadrature
  ξs::Vector{SVector{D, Ftype}}
  ws::Vector{Ftype}
  # shape functions
  Ns::Vector{SVector{N, Ftype}}
  ∇N_ξs::Vector{SMatrix{N, D, Ftype, L}}
end
quadrature_point(e::ReferenceFE, q::Integer) = getfield(e, :ξs)[q]
quadrature_points(e::ReferenceFE) = getfield(e, :ξs)
quadrature_weight(e::ReferenceFE, q::Integer) = getfield(e, :ws)[q]
quadrature_weights(e::ReferenceFE) = getfield(e, :ws)
shape_function_values(e::ReferenceFE) = getfield(e, :Ns)
shape_function_values(e::ReferenceFE, i::Integer) = getfield(e, :Ns)[i]
shape_function_gradients(e::ReferenceFE) = getfield(e, :∇N_ξs)
shape_function_gradients(e::ReferenceFE, i::Integer) = getfield(e, :∇N_ξs)[i]

# shape_function_values(e::ReferenceFE{N, D, Itype, Ftype}, i::Integer) where {N, D, Itype, Ftype} = getindex(getfield(e, :Ns), i)
vertex_nodes(::ReferenceFE{N, D, L, Itype, Ftype}) where {N, D, L, Itype, Ftype} = 1:N

function ReferenceFE(
  e::ReferenceFEType{N, D}, degree::Integer,
  ::Type{Itype} = Int64, ::Type{Ftype} = Float64
) where {N, D, Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates, face_nodes, interior_nodes = element_stencil(e, degree, Itype, Ftype)
  ξs, ws = quadrature_points_and_weights(e, degree, Ftype)
  Ns = Vector{SVector{N, Ftype}}(undef, length(ξs))
  ∇N_ξs = Vector{SMatrix{N, D, Ftype, N * D}}(undef, length(ξs))
  for (q, ξ) in enumerate(ξs)
    Ns[q] = shape_function_values(e, ξ)
    ∇N_ξs[q] = shape_function_gradients(e, ξ)
  end
  return ReferenceFE{N, D, N * D, Itype, Ftype}(
    nodal_coordinates, face_nodes, interior_nodes,
    ξs, ws,
    Ns, ∇N_ξs
  )
end

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
# include("implementations/Tet10.jl")
include("implementations/Tri3.jl")
include("implementations/Tri6.jl")

# include("implementations/SimplexTri.jl")

# precompilation
@setup_workload begin
  @compile_workload begin
    # methods to precompile for all elements
    for el_type in subtypes(ReferenceFEType)
      for degree in [1, 2]
        ReferenceFE(el_type(), degree)
      end
    end
  end
end

end # module ReferenceFEs
