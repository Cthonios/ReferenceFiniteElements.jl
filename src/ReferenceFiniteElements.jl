module ReferenceFiniteElements

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
export EnrichedLagrange
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
export num_edges_per_cell
export num_faces_per_cell
export num_vertices_per_cell
export vertex_coordinates

import PrecompileTools: @compile_workload, @setup_workload
using DocStringExtensions
using FastGaussQuadrature
using LinearAlgebra
using Polynomials
using SpecialPolynomials
using StaticArrays

# includes
include("PolynomialTypes.jl")
include("AbstractElements.jl")
include("QuadratureRules.jl")
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
include("elements/Enriched.jl")

# just precompiling some of the more standard elements for now
# to reduce precompile time...
@setup_workload begin
    @compile_workload begin
        # quadrature rules
        for q in 1:10
            q_rule = GaussLegendre(q)
            q_rule = GaussLobattoLegendre(q)
        end

        # edge elements
        ref_fe = ReferenceFE(Edge{Lagrange, 1}(), GaussLegendre(2))
        ref_fe = ReferenceFE(Edge{Lagrange, 2}(), GaussLegendre(2))

        ref_fe = ReferenceFE(Edge{Lagrange, 1}(), GaussLobattoLegendre(2))
        ref_fe = ReferenceFE(Edge{Lagrange, 2}(), GaussLobattoLegendre(2))

        # hex elements
        ref_fe = ReferenceFE(Hex{Lagrange, 1}(), GaussLegendre(2))
        ref_fe = ReferenceFE(Hex{Lagrange, 1}(), GaussLobattoLegendre(2))

        # quad elements
        ref_fe = ReferenceFE(Quad{Lagrange, 1}(), GaussLegendre(2))
        ref_fe = ReferenceFE(Quad{Lagrange, 2}(), GaussLegendre(2))
        ref_fe = ReferenceFE(Quad{Lagrange, 1}(), GaussLobattoLegendre(2))
        ref_fe = ReferenceFE(Quad{Lagrange, 2}(), GaussLobattoLegendre(2))

        # tet elements
        ref_fe = ReferenceFE(Tet{Lagrange, 1}(), GaussLegendre(2))
        ref_fe = ReferenceFE(Tet{Lagrange, 2}(), GaussLegendre(2))
        ref_fe = ReferenceFE(Tet{Lagrange, 1}(), GaussLobattoLegendre(2))
        ref_fe = ReferenceFE(Tet{Lagrange, 2}(), GaussLobattoLegendre(2))

        ref_fe = ReferenceFE(Tet{EnrichedLagrange, 2}(), GaussLegendre(5))
        ref_fe = ReferenceFE(Tet{EnrichedLagrange, 2}(), GaussLobattoLegendre(2))

        # tri elements
        ref_fe = ReferenceFE(Tri{Lagrange, 1}(), GaussLegendre(2))
        ref_fe = ReferenceFE(Tri{Lagrange, 2}(), GaussLegendre(2))
        ref_fe = ReferenceFE(Tri{Lagrange, 1}(), GaussLobattoLegendre(2))
        ref_fe = ReferenceFE(Tri{Lagrange, 2}(), GaussLobattoLegendre(2))

        ref_fe = ReferenceFE(Tri{EnrichedLagrange, 2}(), GaussLegendre(3))
    end
end

end # module ReferenceFiniteElements
