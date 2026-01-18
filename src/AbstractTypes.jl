abstract type AbstractPolynomialType end

# TODO move somewhere else eventually
struct Hermite <: AbstractPolynomialType
end
struct Lagrange <: AbstractPolynomialType
end
# special type for Vertex
struct NoInterpolation <: AbstractPolynomialType
end
struct RaviartThomas <: AbstractPolynomialType
end
struct Serendipity <: AbstractPolynomialType
end

abstract type AbstractQuadratureType end

cell_quadrature_degree(q_rule::AbstractQuadratureType) = q_rule.cell_degree
surface_quadrature_degree(q_rule::AbstractQuadratureType) = q_rule.surf_degree

struct GaussLobattoLegendre <: AbstractQuadratureType
  cell_degree::Int
  surf_degree::Int

  function GaussLobattoLegendre(degree::Int)
    @assert degree > 0 "Quadrature degree must be greater than zero"
    new(degree, degree)
  end
  
  function GaussLobattoLegendre(cell_degree::Int, surf_degree::Int)
    @assert cell_degree > 0 "Cell quadrature degree must be greater than zero"
    @assert surf_degree > 0 "Surface quadrature degree must be greater than zero"
    new(cell_degree, surf_degree)
  end
end

# methods to define for quadrature types
function cell_quadrature_points_and_weights end
function surface_quadrature_points_and_weights end

"""
$(TYPEDEF)
Base type for the package that all element types
are subtyped off of. The parameters have the following
general meaning.

``PT`` - Polynomial type

``PD`` - Cell polynomial degree

"""
abstract type AbstractElementType{
  PT <: AbstractPolynomialType, 
  PD
} end

"""
$(TYPEDSIGNATURES)
Returns the polynomial type ``PT``.
"""
polynomial_type(::AbstractElementType{PT, PD}) where {PT, PD} = PT
"""
$(TYPEDSIGNATURES)
Returns the polynomial degree ``CD``.
"""
polynomial_degree(::AbstractElementType{PT, PD}) where {PT, PD} = PD

# methods to define for elements for compile time info

# topology methods
function boundary_element end
function boundary_normals end # always returns 3 x num boundaries
function dimension end
function edge_vertices end
function face_vertices end
function num_boundaries end
function num_edges end
function num_faces end
function num_vertices_per_cell end
# function num_vertices_per_edge end
# function num_vertices_per_face end
function vertex_coordinates end # always returns 3 x num_vertices_per_cell

# derived topology methods
cell_vertices(e::AbstractElementType) = 1:num_vertices_per_cell(e) |> collect

# dof related methods
# TODO should seperate dofs into
# 1. vertex dofs
# 2. edge dofs
# 3. face dofs
# 4. cell dofs (what we're calling interior dofs right now)
# this is how defelement defines things
function boundary_dofs end
function dof_coordinates end
function interior_dofs end
function num_cell_dofs end
# below method does not include e.g. vertices
function num_dofs_on_boundary end
function num_interior_dofs end


# 0d elements
"""
$(TYPEDEF)
"""
abstract type AbstractVertex <: AbstractElementType{NoInterpolation, 0} end

boundary_element(::AbstractVertex, ::Int) = nothing
boundary_normals(::AbstractVertex) = Matrix{Float64}(undef, 3, 0)
dimension(::AbstractVertex) = 0
edge_vertices(::AbstractVertex) = Matrix{Int}(undef, 0, 0)
face_vertices(::AbstractVertex) = Matrix{Int}(undef, 0, 0)
num_boundaries(::AbstractVertex) = 0
num_edges(::AbstractVertex) = 0
num_faces(::AbstractVertex) = 0
num_vertices_per_cell(::AbstractVertex) = 1
# num_vertices_per_edge(e::AbstractVertex) = 0
# num_vertices_per_face(::AbstractVertex) = 0
vertex_coordinates(::AbstractVertex) = [0. 0. 0.]' |> collect

# this one is an oddity
boundary_dofs(::AbstractVertex) = Matrix{Int}(undef, 0, 0)
dof_coordinates(::AbstractVertex) = [0. 0. 0.]' |> collect
interior_dofs(::AbstractVertex) = Vector{Int}(undef, 0)
num_cell_dofs(::AbstractVertex) = 1
num_dofs_on_boundary(::AbstractVertex, ::Int) = 0
num_interior_dofs(::AbstractVertex) = 0

# 1d elements
"""
$(TYPEDEF)
"""
abstract type AbstractEdge{PT, PD} <: AbstractElementType{PT, PD} end
boundary_element(::AbstractEdge, ::Int) = Vertex()
boundary_normals(::AbstractEdge) = [-1. 0. 0.; 1. 0. 0.]' |> collect
dimension(::AbstractEdge) = 1
edge_vertices(::AbstractEdge) = [1 2]' |> collect
face_vertices(::AbstractEdge) = Matrix{Int}(undef, 0, 0)
num_boundaries(::AbstractEdge) = 2
num_edges(::AbstractEdge) = 1
num_faces(::AbstractEdge) = 0
num_vertices_per_cell(::AbstractEdge) = 2
# num_vertices_per_edge(::AbstractEdge) = 2
# num_vertices_per_face(::AbstractEdge) = 0
function vertex_coordinates(e::AbstractEdge) 
  if e.shifted
    return [
      0. 1.;
      0. 0.;
      0. 0.
    ]
  else
    return [
      -1. 1.;
       0. 0.;
       0. 0.
    ]
  end
end

# 2d elements
"""
$(TYPEDEF)
"""
abstract type AbstractFace{PT, PD} <: AbstractElementType{PT, PD} end
dimension(::AbstractFace) = 2
face_vertices(e::AbstractFace) = 1:num_vertices_per_cell(e) |> collect
num_boundaries(e::AbstractFace) = num_edges(e)
num_faces(::AbstractFace) = 1
# num_vertices_per_edge(e::AbstractFace) = 2
# num_vertices_per_face(e::AbstractFace) = num_vertices_per_cell(e)

"""
$(TYPEDEF)
"""
abstract type AbstractQuad{PT, PD} <: AbstractFace{PT, PD} end
boundary_element(::AbstractQuad{PT, PD}, ::Int) where {PT, PD} = Edge{PT, PD}()
boundary_normals(::AbstractQuad) = [
  0. 1. 0. -1.;
 -1. 0. 1.  0.;
  0. 0. 0.  0.
]
edge_vertices(::AbstractQuad) = [
  1 2 3 4;
  2 3 4 1
]
num_edges(::AbstractQuad) = 4
num_vertices_per_cell(::AbstractQuad) = 4
vertex_coordinates(::AbstractQuad) = [
  -1.  1. 1. -1.;
  -1. -1. 1.  1.;
   0.  0. 0.  0.
]

"""
$(TYPEDEF)
"""
abstract type AbstractTri{PT, PD} <: AbstractFace{PT, PD} end
boundary_element(::AbstractTri{PT, PD}, ::Int) where {PT, PD} = Edge{PT, PD}(; shifted = true)
boundary_normals(::AbstractTri) = [
  0. 1. / sqrt(2.) -1.;
 -1. 1. / sqrt(2.)  0.;
  0. 0.             0.
]
edge_vertices(::AbstractTri) = [
  1 2 3;
  2 3 1
]
num_edges(::AbstractTri) = 3
num_vertices_per_cell(::AbstractTri) = 3
vertex_coordinates(::AbstractTri) = [
  0. 1. 0.;
  0. 0. 1.;
  0. 0. 0.
]

# TODO finish this abstract interface
# 3d elements
"""
$(TYPEDEF)
"""
abstract type AbstractVolume{PT, PD} <: AbstractElementType{PT, PD} end
dimension(::AbstractVolume) = 3
num_boundaries(e::AbstractVolume) = num_faces(e)
# num_vertices_per_edge(e::AbstractVolume) = 2
# num_vertices_per_face(e::AbstractVolume) = num_vertices_per_cell(boundary_element(e))

"""
$(TYPEDEF)
"""
abstract type AbstractHex{PT, PD} <: AbstractVolume{PT, PD} end
boundary_element(::AbstractHex{PT, PD}, ::Int) where {PT, PD} = Quad{PT, PD}()
boundary_normals(::AbstractHex) = [
   0. 1. 0. -1.  0. 0.
   0. 0. 0.  0. -1. 1.
  -1. 0. 1.  0.  0. 0.
]
edge_vertices(::AbstractHex) = [
  1 2 3 4 5 6 7 8 1 2 3 4
  2 3 4 1 6 7 8 5 5 6 7 8
]
face_vertices(::AbstractHex) = [
  1 2 3 1 1 5
  2 3 4 5 4 6
  6 7 8 8 3 7
  5 6 7 4 2 8
]
num_edges(::AbstractHex) = 12
num_faces(::AbstractHex) = 6
num_vertices_per_cell(::AbstractHex) = 8
vertex_coordinates(::AbstractHex) = [
  -1.  1.  1. -1. -1.  1. 1. -1.
  -1. -1.  1.  1. -1. -1. 1.  1.
  -1. -1. -1. -1.  1.  1. 1.  1.
]

"""
$(TYPEDEF)
"""
abstract type AbstractPyramid{PT, PD} <: AbstractVolume{PT, PD} end
# TODO finish me

"""
$(TYPEDEF)
"""
abstract type AbstractTet{PT, PD} <: AbstractVolume{PT, PD} end
boundary_element(::AbstractTet{PT, PD}, ::Int) where {PT, PD} = Tri{PT, PD}()
boundary_normals(::AbstractTet) = [
   0. 1. / sqrt(3.) -1.  0.
  -1. 1. / sqrt(3.)  0.  0.
   0. 1. / sqrt(3.)  0. -1.
]
edge_vertices(::AbstractTet) = [
  1 2 3 1 2 3
  2 3 1 4 4 4
]
face_vertices(::AbstractTet) = [
  1 2 1 1
  2 3 4 3
  4 4 3 2
]
num_edges(::AbstractTet) = 6
num_faces(::AbstractTet) = 4
num_vertices_per_cell(::AbstractTet) = 4
vertex_coordinates(::AbstractTet) = [
  0. 1. 0. 0.
  0. 0. 1. 0.
  0. 0. 0. 1.
]

"""
$(TYPEDEF)
"""
abstract type AbstractWedge{PT, PD} <: AbstractVolume{PT, PD} end
# TODO finish out