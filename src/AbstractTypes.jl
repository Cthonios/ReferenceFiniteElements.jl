"""
$(TYPEDEF)
Base type for the package that all element types
are subtyped off of. The parameters have the following
general meaning.

``D`` - Dimension

``V`` - Number of vertices

``E`` - Number of edges

``F`` - Number of faces

``I`` - Interpolation type

``P`` - Polynomial degree

``Q`` - quadrature degree

"""
abstract type AbstractElementType{D, V, E, F, I, P, Q} end

"""
$(TYPEDSIGNATURES)
Returns dimension D.
"""
dimension(::AbstractElementType{D, V, E, F, I, P, Q}) where {D, V, E, F, I, P, Q} = D
"""
$(TYPEDSIGNATURES)
"""
num_vertices(::AbstractElementType{D, V, E, F, I, P, Q}) where {D, V, E, F, I, P, Q} = V
"""
$(TYPEDSIGNATURES)
"""
num_edges(::AbstractElementType{D, V, E, F, I, P, Q}) where {D, V, E, F, I, P, Q} = E
"""
$(TYPEDSIGNATURES)
"""
num_faces(::AbstractElementType{D, V, E, F, I, P, Q}) where {D, V, E, F, I, P, Q} = F
"""
$(TYPEDSIGNATURES)
Returns the interpolation type ``I``.
"""
interpolation_type(::AbstractElementType{D, V, E, F, I, P, Q}) where {D, V, E, F, I, P, Q} = I
"""
$(TYPEDSIGNATURES)
Returns the interpolation polynomial degree ``P``.
"""
polynomial_degree(::AbstractElementType{D, V, E, F, I, P, Q}) where {D, V, E, F, I, P, Q} = P
"""
$(TYPEDSIGNATURES)
Returns the interpolation polynomial degree ``P``.
"""
polynomial_degree(::Type{<:AbstractElementType{D, V, E, F, I, P, Q}}) where {D, V, E, F, I, P, Q} = P
"""
$(TYPEDSIGNATURES)
Return the quadrature degree ``Q``
"""
quadrature_degree(::AbstractElementType{D, V, E, F, I, P, Q}) where {D, V, E, F, I, P, Q} = Q
"""
$(TYPEDSIGNATURES)
"""
function num_interior_vertices end
"""
$(TYPEDSIGNATURES)
"""
function num_quadrature_points end
"""
$(TYPEDSIGNATURES)
"""
function num_vertices_per_edge end
"""
$(TYPEDSIGNATURES)
"""
function num_vertices_per_face end
# """
# $(TYPEDSIGNATURES)
# Returns the number of edges in the element topology.
# """
# num_edges(e::AbstractElementType) = num_edges(typeof(e))
# """
# $(TYPEDSIGNATURES)
# Returns the number of faces in the element topology.
# """
# num_faces(e::AbstractElementType) = num_faces(typeof(e))
# """
# $(TYPEDSIGNATURES)
# Returns the number of interior vertices in the element topology.
# """
# num_interior_vertices(e::AbstractElementType) = num_interior_vertices(typeof(e))
# """
# $(TYPEDSIGNATURES)
# Returns the number of quadrature points.
# """
# num_quadrature_points(e::AbstractElementType) = num_quadrature_points(typeof(e))
# """
# $(TYPEDSIGNATURES)
# Returns the number of vertices in the element topology.
# """
# num_vertices(e::AbstractElementType) = num_vertices(typeof(e))
# """
# $(TYPEDSIGNATURES)
# Returns the number of vertices per edge in the element topology.
# TODO this may not generalize well to arbitrary elements
# """
# num_vertices_per_edge(e::AbstractElementType) = num_vertices_per_edge(e)
# """
# $(TYPEDSIGNATURES)
# Returns the number of vertices per face in the element topology.
# TODO this may not generalize well to arbitrary elements
# """
# num_vertices_per_face(e::AbstractElementType) = num_vertices_per_face(e)

# 0d elements
"""
$(TYPEDEF)
"""
abstract type AbstractVertex{I, P, Q} <: AbstractElementType{0, 1, 0, 0, I, P, Q} end
num_interior_vertices(::AbstractVertex) = 0
num_quadrature_points(::AbstractVertex) = 1
num_vertices_per_edge(::AbstractVertex) = 0
num_vertices_per_face(e::AbstractVertex) = 0

# 1d elements
"""
$(TYPEDEF)
"""
abstract type AbstractEdge{V, I, P, Q} <: AbstractElementType{1, V, 1, 0, I, P, Q} end
function num_interior_vertices(e::AbstractEdge{V, I, P, Q}) where {V, I, P, Q}
  if P < 2
    return 0
  else
    return polynomial_degree(e) - 1
  end
end
num_quadrature_points(::AbstractEdge{V, I, P, Q}) where {V, I, P, Q} = Q
num_vertices_per_edge(e::AbstractEdge) = num_vertices(e)
num_vertices_per_face(e::AbstractEdge) = 0

# 2d elements
"""
$(TYPEDEF)
"""
abstract type AbstractFace{V, E, I, P, Q} <: AbstractElementType{2, V, E, 1, I, P, Q} end
num_vertices_per_edge(e::AbstractFace) = num_vertices(surface_element(e))
num_vertices_per_face(e::AbstractFace) = num_vertices(e)

"""
$(TYPEDEF)
"""
abstract type AbstractQuad{V, I, P, Q} <: AbstractFace{V, 4, I, P, Q} end
num_interior_vertices(e::AbstractQuad) = num_interior_vertices(surface_element(e))^2
num_quadrature_points(e::AbstractQuad) = num_quadrature_points(surface_element(e))^2

"""
$(TYPEDEF)
"""
abstract type AbstractTri{V, I, P, Q} <: AbstractFace{V, 3, I, P, Q} end
num_interior_vertices(e::AbstractTri) = num_interior_vertices(surface_element(e)) * (num_interior_vertices(surface_element(e)) + 1) ÷ 2
num_quadrature_points(e::AbstractTri) = num_quadrature_points(surface_element(e)) * (num_quadrature_points(surface_element(e)) + 1) ÷ 2

# 3d elements
"""
$(TYPEDEF)
"""
abstract type AbstractVolume{V, E, F, I, P, Q} <: AbstractElementType{3, V, E, F, I, P, Q} end
num_vertices_per_edge(e::AbstractVolume) = num_vertices(surface_element(surface_element(e)))
num_vertices_per_face(e::AbstractVolume) = num_vertices(surface_element(e))

"""
$(TYPEDEF)
"""
abstract type AbstractHex{V, I, P, Q} <: AbstractVolume{V, 12, 6, I, P, Q} end
num_interior_vertices(e::AbstractHex) = num_interior_vertices(surface_element(surface_element(e)))^3
num_quadrature_points(e::AbstractHex) = num_quadrature_points(surface_element(surface_element(e)))^3

"""
$(TYPEDEF)
"""
abstract type AbstractTet{V, I, P, Q} <: AbstractVolume{V, 6, 4, I, P, Q} end
function num_interior_vertices(e::AbstractTet{V, I, P, Q}) where {V, I, P, Q}
  if P < 3
    return 0
  else
    # TODO check this
    return num_vertices(e) - num_faces(e) * num_vertices_per_edge(e) + 4
  end
end
function num_quadrature_points(e::AbstractTet)
  n_q_edge = num_quadrature_points(surface_element(surface_element(e)))
  return (n_q_edge + 1) * (n_q_edge + 2) * (n_q_edge + 3) ÷ 6
end

# TODO add pyramid and wedge type

# containers for interpolants
"""
$(TYPEDEF)
"""
abstract type AbstractInterpolationType end
"""
$(TYPEDEF)
"""
abstract type AbstractInterpolantsContainer{I} end

# interpolants stuff
"""
$(TYPEDSIGNATURES)
"""
function num_shape_functions(e::AbstractElementType{D, V, E, F, Lagrange, P, Q}) where {D, V, E, F, P, Q}
  if P == 0
    return 1
  else
    return num_vertices(e)
  end
end
