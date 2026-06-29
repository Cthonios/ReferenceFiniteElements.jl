"""
$(TYPEDEF)
"""
abstract type AbstractPolynomialType end
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
