module ReferenceFiniteElements

# import Polynomials: derivative
# import SpecialPolynomials: basis, Legendre, ShiftedLegendre
using DocStringExtensions
using FastGaussQuadrature
using LinearAlgebra
using Polynomials
using SpecialPolynomials
using StaticArrays

function Polynomials.derivative(p::P) where {
    B <: SpecialPolynomials.ShiftedLegendreBasis,
    T, X,
    P <: SpecialPolynomials.AbstractUnivariatePolynomial{B,T,X}
}
    hasnan(p) && return Polynomials.⟒(P){T,X}(T[NaN])

    d = degree(p)
    qs = zeros(T, d)

    for i in 0:(d - 1)
        gamma = 2 * (2i + 1)   # ← extra factor of 2 from chain rule
        qs[i + 1] = gamma * sum(p[j] for j in (i + 1):2:d)
    end

    dp = Polynomials.MutableDensePolynomial{B, Float64, :x}(qs)
    return dp
end

# includes
include("AbstractTypes.jl")
# include("ArrayBackends.jl")
# include("Interpolants.jl")
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
include("elements/Tet.jl")

# abstract interface
export polynomial_degree
export polynomial_type

# dof interface
export boundary_dofs
export dof_coordinates
export interior_dofs
export num_cell_dofs
export num_interior_dofs
# export surface_dof_coordinates

# elements
export Edge
export Hex
export Quad
export Tet
export Tri
export Vertex
# export Edge0, Edge2, Edge3, Edge
# export Hex0, Hex8, Hex
# export Quad0, Quad4, Quad9, Quad
# export Tet0, Tet4, Tet10
# export Tri0, Tri3, Tri6, Tri
# export Vertex

# polynomial types
export Lagrange

# quadrature types
export GaussLobattoLegendre

# reference fe
export ReferenceFE
export cell_quadrature_point
export cell_quadrature_weight
export cell_shape_function_gradient
export cell_shape_function_hessian
export cell_shape_function_value
export num_cell_quadrature_points
export num_surface_quadrature_points
export surface_quadrature_point
export surface_quadrature_weight
export surface_shape_function_gradient
export surface_shape_function_hessian
export surface_shape_function_value

# topology interface
export boundary_element
export boundary_normals
export cell_vertices
export dimension
export edge_vertices
export face_vertices
export num_boundaries
export num_edges
export num_faces
export num_vertices_per_cell
export num_vertices_per_edge
export num_vertices_per_face
export vertex_coordinates

# # main type
# export ReferenceFE,
#        MappedInterpolants,
#        MappedSurfaceInterpolants

# # methods
# export num_faces,
#        num_quadrature_points,
#        num_shape_functions,
#        num_vertices,
#        surface_element

# export quadrature_point,
#        quadrature_points,
#        quadrature_weight,
#        quadrature_weights,
#        shape_function_gradient,
#        shape_function_gradients,
#        shape_function_hessian,
#        shape_function_hessians,
#        shape_function_value,
#        shape_function_values,
#        surface_normal,
#        surface_normals,
#        surface_quadrature_point,
#        surface_quadrature_points,
#        surface_quadrature_weight,
#        surface_quadrature_weights,
#        surface_shape_function_gradient,
#        surface_shape_function_gradients,
#        surface_shape_function_hessian,
#        surface_shape_function_hessians,
#        surface_shape_function_value,
#        surface_shape_function_values

# export calculate_JxW,
#        integrate,
#        map_shape_function_gradient,
#        mapping_jacobian

end # module ReferenceFiniteElements
