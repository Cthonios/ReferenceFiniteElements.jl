module ReferenceFiniteElements

# element types
# export Edge
export Hex8,
       Quad4, Quad9,
       SimplexTri,
       Tet4, Tet10,
       Tri3, Tri6

# types
export ReferenceFE

# methods.
export element_type,
       num_dimensions,
       num_nodes_per_element,
       num_q_points,
       quadrature_points,
       quadrature_weights,
       shape_function_gradients,
       shape_function_hessians,
       shape_function_values,
       surface_quadrature_points,
       surface_quadrature_weights,
       surface_shape_function_gradients,
       surface_shape_function_hessians,
       surface_shape_function_values,
       vertex_nodes

# dependencies
using DocStringExtensions,
      FastGaussQuadrature,
      LinearAlgebra,
      Polynomials,
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

# Types
include("ReferenceFETypes.jl")
include("Interpolants.jl")
include("ReferenceFEs.jl")
# include("Utils.jl")

# implementations of things common across multiple element types
include("implementations/HexCommon.jl")
include("implementations/QuadCommon.jl")
include("implementations/TetCommon.jl")
include("implementations/TriCommon.jl")

# implementations of things specific to element types
include("implementations/Edge.jl")
include("implementations/Hex8.jl")
include("implementations/Quad4.jl")
include("implementations/Quad9.jl")
include("implementations/SimplexTri.jl")
include("implementations/Tet4.jl")
include("implementations/Tet10.jl")
include("implementations/Tri3.jl")
include("implementations/Tri6.jl")

end # module ReferenceFEs
