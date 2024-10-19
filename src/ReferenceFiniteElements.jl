module ReferenceFiniteElements

using DocStringExtensions
using FastGaussQuadrature
using LinearAlgebra
using Polynomials
using SpecialPolynomials
using StaticArrays
using StructArrays

include("AbstractTypes.jl")
include("ArrayBackends.jl")
include("Interpolants.jl")
include("ReferenceFEs.jl")

# 0-d elements
include("elements/Vertex.jl")

# 1-d elements
include("elements/Edge.jl")

# 2-d elements
include("elements/Quad.jl")
include("elements/Tri.jl")

# 3-d elements
include("elements/Hex.jl")
# include("elements/Tet.jl")

# elements
export Edge0, Edge2, Edge3, Edge
export Hex0, Hex8, Hex
export Quad0, Quad4, Quad9, Quad
export Tri0, Tri3, Tri6, Tri
export Vertex

# interpolations
export Lagrange

# main type
export ReferenceFE

# methods
export num_faces,
       num_nodes,
       num_quadrature_points,
       num_shape_functions,
       surface_element

export quadrature_point,
       quadrature_points,
       quadrature_weight,
       quadrature_weights,
       shape_function_gradient,
       shape_function_gradients,
       shape_function_hessian,
       shape_function_hessians,
       shape_function_value,
       shape_function_values,
       surface_normal,
       surface_normals,
       surface_quadrature_point,
       surface_quadrature_points,
       surface_quadrature_weight,
       surface_quadrature_weights,
       surface_shape_function_gradient,
       surface_shape_function_gradients,
       surface_shape_function_hessian,
       surface_shape_function_hessians,
       surface_shape_function_value,
       surface_shape_function_values

export calculate_JxW,
       integrate,
       map_shape_function_gradient,
       mapping_jacobian

end # module ReferenceFiniteElements
