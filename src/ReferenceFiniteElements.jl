module ReferenceFiniteElements

# element types
export Edge
export Hex8
export Quad4, Quad9
export Tet4, Tet10
export Tri3, Tri6

# types
export Quadrature
export ReferenceFE
export ReferenceFEStencil
export ShapeFunctions

# methods
export num_dimensions
export num_q_points
export quadrature_points
export quadrature_weights
export reference_fe_coordinates
export reference_fe_face_nodes
export reference_fe_interior_nodes
export reference_fe_type
export reference_fe_vertex_nodes
export shape_function_gradient
export shape_function_gradients
export shape_function_values

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

# interface
include("ElementStencil.jl")
include("Quadrature.jl")
include("ShapeFunctions.jl")

struct ReferenceFE{N, D, Itype <: Integer, Rtype <: Real, RefFE <: ReferenceFEType}
  q_rule::Quadrature{Rtype}
  stencil::ReferenceFEStencil{Itype, Rtype, RefFE}
  shape_functions::ShapeFunctions{N, D, Rtype}
end

function ReferenceFE(
  e::ReferenceFEType{N, D}, degree::Integer, 
  Itype = Integer, Rtype = Float64
) where {N, D}
  q_rule = Quadrature(e, degree)
  stencil = ReferenceFEStencil(e, degree)
  shape_functions = ShapeFunctions(e, q_rule)
  return ReferenceFE{N, D, Itype, Rtype, ReferenceFEType{N, D}}(q_rule, stencil, shape_functions)
end

# implementations
include("implementations/HexCommon.jl")
include("implementations/QuadCommon.jl")
include("implementations/TetCommon.jl")
include("implementations/TriCommon.jl")

include("implementations/Edge.jl")
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
    # methods to precompile for all elements
    for el_type in subtypes(ReferenceFEType)
      for degree in [1, 2]
        ReferenceFE(el_type(), degree)
      end
      # for el_type in subtypes(abstract_el_type)
      #   for degree in [1, 2]
      #     ReferenceFE(el_type(), degree)
      #   end
      # end
    end
  end
end

end # module ReferenceFEs
