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
# export num_dimensions
# export num_q_points
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
struct ReferenceFE{N, D, Ndof, Itype <: Integer, Ftype <: AbstractFloat}
  # element nodal info
  nodal_coordinates::VecOrMat{Ftype}
  face_nodes::Matrix{Itype}
  interior_nodes::Vector{Itype}
  # quadrature
  ξs::VecOrMat{Ftype}
  ws::Vector{Ftype}
  # shape functions
  Ns::Vector{SVector{N, Ftype}}
  ∇N_ξs::Vector{SMatrix{N, D, Ftype, Ndof}}
end
vertex_nodes(::ReferenceFE{N, D, Itype, Ftype}) where {N, D, Itype, Ftype} = 1:N

function ReferenceFE(
  e::ReferenceFEType{N, D}, degree::Integer,
  ::Type{Itype} = Int64, ::Type{Ftype} = Float64
) where {N, D, Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates, face_nodes, interior_nodes = element_stencil(e, degree, Itype, Ftype)
  ξs, ws = quadrature_points_and_weights(e, degree, Ftype)
  Ns = shape_function_values.((e,), eachcol(ξs))
  ∇N_ξs = shape_function_gradients.((e,), eachcol(ξs))
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
