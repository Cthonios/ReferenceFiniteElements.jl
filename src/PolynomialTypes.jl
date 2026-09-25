"""
$(TYPEDEF)
"""
abstract type AbstractPolynomialType end
"""
$(TYPEDEF)
Lagrange functions of degree ``PD`` enriched with bubble functions, carried by
a nodal basis (unit value at one node, zero at the others). Implemented at
degree 2 for the triangle, the six ``P_2`` functions and the cubic interior
bubble on seven nodes, and for the tetrahedron, the ten ``P_2`` functions, one
cubic bubble per face and the quartic interior bubble on fifteen nodes, which
is the conforming three-dimensional Crouzeix--Raviart displacement space.
"""
struct EnrichedLagrange <: AbstractPolynomialType
end
"""
$(TYPEDEF)
"""
struct Hermite <: AbstractPolynomialType
end
"""
$(TYPEDEF)
"""
struct Lagrange <: AbstractPolynomialType
end
"""
$(TYPEDEF)
Planned
"""
struct NedelecFirstKind <: AbstractPolynomialType
end
"""
$(TYPEDEF)
Planned
"""
struct NedelecSecondKind <: AbstractPolynomialType
end
"""
$(TYPEDEF)
special type for ``Vertex``
"""
struct NoInterpolation <: AbstractPolynomialType
end
"""
$(TYPEDEF)
"""
struct RaviartThomas <: AbstractPolynomialType
end
"""
$(TYPEDEF)
"""
struct Serendipity <: AbstractPolynomialType
end
