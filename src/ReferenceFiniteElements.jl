module ReferenceFiniteElements

using Adapt
using DocStringExtensions
using FastGaussQuadrature
using LinearAlgebra
using Polynomials
using SpecialPolynomials
using StaticArrays

# includes
include("AbstractTypes.jl")
include("ReferenceFE.jl")

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
export num_dofs_on_boundary
export num_interior_dofs
# export surface_dof_coordinates

# elements
export Edge
export Hex
export Quad
export Tet
export Tri
export Vertex

# interpolants
# export H1OrL2Interpolants
export MappedH1OrL2Interpolants
export MappedH1OrL2SurfaceInterpolants
export StaticH1OrL2Interpolants
export StaticH1OrL2InterpolantsWithHessians

# polynomial types
export Hermite
export Lagrange
export RaviartThomas
export Serendipity

# quadrature types
export AbstractQuadratureType
export GaussLegendre
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
export boundary_normal
export boundary_normals
export cell_vertices
export dimension
export edge_vertices
export element
export face_vertices
export num_boundaries
export num_edges
export num_faces
export num_vertices_per_cell
export vertex_coordinates

end # module ReferenceFiniteElements
