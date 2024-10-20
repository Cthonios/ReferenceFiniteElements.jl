# D - dimension
# I - Interpolation type
# P - polynomial
# R - quadrature rule
# Q - quadrature order
abstract type AbstractElementType{D, I, P, Q} end

dimension(::AbstractElementType{D, I, P, Q}) where {D, I, P, Q} = D
interpolation_type(::AbstractElementType{D, I, P, Q}) where {D, I, P, Q} = I
polynomial_degree(::AbstractElementType{D, I, P, Q}) where {D, I, P, Q} = P
polynomial_degree(::Type{<:AbstractElementType{D, I, P, Q}}) where {D, I, P, Q} = P
quadrature_degree(::AbstractElementType{D, I, P, Q}) where {D, I, P, Q} = Q

num_edges(e::AbstractElementType) = num_edges(typeof(e))
num_faces(e::AbstractElementType) = num_faces(typeof(e))
num_interior_vertices(e::AbstractElementType) = num_interior_vertices(typeof(e))
num_quadrature_points(e::AbstractElementType) = num_quadrature_points(typeof(e))
num_vertices(e::AbstractElementType) = num_vertices(typeof(e))
num_vertices_per_edge(e::AbstractElementType) = num_vertices_per_edge(typeof(e))
num_vertices_per_face(e::AbstractElementType) = num_vertices_per_face(typeof(e))

# 0d elements
abstract type AbstractVertex{I, P, Q} <: AbstractElementType{0, I, P, Q} end
num_edges(::Type{<:AbstractVertex}) = 0
num_faces(::Type{<:AbstractVertex}) = 0
num_interior_vertices(::Type{<:AbstractVertex}) = 0
num_quadrature_points(::Type{<:AbstractVertex}) = 1
num_vertices(::Type{<:AbstractVertex}) = 1
num_vertices_per_edge(::Type{<:AbstractVertex}) = 0
num_vertices_per_face(e::Type{<:AbstractVertex}) = 0
surface_element_type(::Type{<:AbstractVertex}) = nothing

# 1d elements
abstract type AbstractEdge{I, P, Q} <: AbstractElementType{1, I, P, Q} end
num_edges(::Type{<:AbstractEdge}) = 1
num_faces(::Type{<:AbstractEdge}) = 0
function num_interior_vertices(e::Type{<:AbstractEdge{I, P, Q}}) where {I, P, Q}
  if P < 2
    return 0
  else
    return polynomial_degree(e) - 1
  end
end
num_quadrature_points(::Type{<:AbstractEdge{I, P, Q}}) where {I, P, Q} = Q
function num_vertices(e::Type{<:AbstractEdge{I, P, Q}}) where {I, P, Q}
  if P == 0
    return 2
  else
    return polynomial_degree(e) + 1
  end
end
num_vertices_per_edge(e::Type{<:AbstractEdge}) = num_vertices(e)
num_vertices_per_face(e::Type{<:AbstractEdge}) = 0
surface_element_type(::Type{<:AbstractEdge}) = AbstractVertex

# 2d elements
abstract type AbstractFace{I, P, Q} <: AbstractElementType{2, I, P, Q} end
num_faces(e::Type{<:AbstractFace}) = 1
num_vertices_per_edge(e::Type{<:AbstractFace}) = num_vertices(surface_element_type(e))
num_vertices_per_face(e::Type{<:AbstractFace}) = num_vertices(e)
surface_element_type(::Type{<:AbstractFace{I, P, Q}}) where {I, P, Q} = AbstractEdge{I, P, Q}

abstract type AbstractQuad{I, P, Q} <: AbstractFace{I, P, Q} end
num_edges(e::Type{<:AbstractQuad}) = 4
num_interior_vertices(e::Type{<:AbstractQuad}) = num_interior_vertices(surface_element_type(e))^2
num_quadrature_points(e::Type{<:AbstractQuad}) = num_quadrature_points(surface_element_type(e))^2
num_vertices(e::Type{<:AbstractQuad}) = num_vertices(surface_element_type(e))^2

abstract type AbstractTri{I, P, Q} <: AbstractFace{I, P, Q} end
num_edges(e::Type{<:AbstractTri}) = 3
num_interior_vertices(e::Type{<:AbstractTri}) = num_interior_vertices(surface_element_type(e)) * (num_interior_vertices(surface_element_type(e)) + 1) ÷ 2
num_quadrature_points(e::Type{<:AbstractTri}) = num_quadrature_points(surface_element_type(e)) * (num_quadrature_points(surface_element_type(e)) + 1) ÷ 2
num_vertices(e::Type{<:AbstractTri}) = num_vertices(surface_element_type(e)) * (num_vertices(surface_element_type(e)) + 1) ÷ 2

# 3d elements
abstract type AbstractVolume{I, P, Q} <: AbstractElementType{3, I, P, Q} end
num_vertices_per_edge(e::Type{<:AbstractVolume}) = num_vertices(surface_element_type(surface_element_type(e)))
num_vertices_per_face(e::Type{<:AbstractVolume}) = num_vertices(surface_element_type(e))

abstract type AbstractHex{I, P, Q} <: AbstractVolume{I, P, Q} end
num_edges(e::Type{<:AbstractHex}) = 12
num_faces(e::Type{<:AbstractHex}) = 6
num_interior_vertices(e::Type{<:AbstractHex}) = num_interior_vertices(surface_element_type(surface_element_type(e)))^3
num_quadrature_points(e::Type{<:AbstractHex}) = num_quadrature_points(surface_element_type(surface_element_type(e)))^3
num_vertices(e::Type{<:AbstractHex}) = num_vertices(surface_element_type(surface_element_type(e)))^3
surface_element_type(::Type{<:AbstractHex{I, P, Q}}) where {I, P, Q} = AbstractQuad{I, P, Q}

abstract type AbstractTet{I, P, Q} <: AbstractVolume{I, P, Q} end
num_edges(e::Type{<:AbstractTet}) = 6
num_faces(e::Type{<:AbstractTet}) = 4
function num_interior_vertices(e::AbstractTet{I, P, Q}) where {I, P, Q}
  if P < 3
    return 0
  else
    # TODO check this
    return num_vertices(e) - num_faces(e) * num_vertices_per_edge(e) + 4
  end
end
function num_quadrature_points(e::Type{<:AbstractTet})
  n_q_edge = num_quadrature_points(surface_element_type(surface_element_type(e)))
  return (n_q_edge + 1) * (n_q_edge + 2) * (n_q_edge + 3) ÷ 6
end
function num_vertices(e::Type{<:AbstractTet})
  if polynomial_degree(e) == 0
    return 4 
  else
    return (polynomial_degree(e) + 1) * (polynomial_degree(e) + 2) * (polynomial_degree(e) + 3) ÷ 6
  end
end
surface_element_type(::Type{<:AbstractTet{I, P, Q}}) where {I, P, Q} = AbstractTri{I, P, Q} 

# TODO add pyramid and wedge type

# containers for interpolants
abstract type AbstractInterpolationType end
abstract type AbstractInterpolantsContainer{I} end

# interpolants stuff
function num_shape_functions(e::AbstractElementType{D, Lagrange, P, Q}) where {D, P, Q}
  if P == 0
    return 1
  else
    return num_vertices(e)
  end
end
