module ReferenceFiniteElements

using FastGaussQuadrature
using InteractiveUtils
using LinearAlgebra
using Polynomials
using PrecompileTools
using SpecialPolynomials
using StaticArrays
using StructArrays

# abstract interface
abstract type AbstractReferenceFE end

# interface
include("ElementStencil.jl")
include("Quadrature.jl")
include("ShapeFunctions.jl")

# function ReferenceFE(e::E, degree::Int) where E <: ReferenceFE
#   q_rule = Quadrature(e, degree)
#   stencil = ElementStencil(e, degree)
#   shape_functions = ShapeFunctions(e, degree)
#   return q_rule, stencil, shape_functions
# end

struct ReferenceFE{N, D, Itype <: Integer, Rtype <: Real, RefFE <: AbstractReferenceFE}
  q_rule::Quadrature{Rtype}
  stencil::ElementStencil{Itype, Rtype, RefFE}
  shape_functions::StructArray{ShapeFunctionPair{N, D, Rtype}}
end

function ReferenceFE(e::E, degree::Integer, Itype = Integer, Rtype = Float64) where E <: AbstractReferenceFE
  q_rule = Quadrature(e, degree, )
  stencil = ElementStencil(e, degree)
  shape_functions = ShapeFunctions(e, degree)
  num_nodes, num_dim = size(shape_functions.∇N_ξ[1])
  return ReferenceFE{num_nodes, num_dim, Itype, Rtype, E}(q_rule, stencil, shape_functions)
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
# @setup_workload begin
#   @compile_workload begin
#     # methods to precompile for all elements
#     for abstract_el_type in subtypes(AbstractReferenceFE)
#       for el_type in subtypes(abstract_el_type)
#         ElementStencil(el_type(), 1)
#         Quadrature(el_type(), 1)
#         ShapeFunctions(el_type(), 1)
#       end
#     end
#   end
# end

export AbstractReferenceFE
export ReferenceFE

end # module ReferenceFEs
