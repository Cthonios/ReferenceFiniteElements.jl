# Abstract types

## Element abstract type hierarchy
The root of the type hierarchy tree in ```ReferenceFiniteElements``` is the following 
abstract type.

```@docs
ReferenceFiniteElements.AbstractElementType
```

## Interpolants abstract types
```@docs
ReferenceFiniteElements.AbstractInterpolants
ReferenceFiniteElements.AbstractPolynomialType
polynomial_degree
polynomial_type
```

## Quadrature abstract types
```@docs
ReferenceFiniteElements.AbstractQuadratureType
ReferenceFiniteElements.cell_quadrature_degree
ReferenceFiniteElements.cell_quadrature_points_and_weights
ReferenceFiniteElements.surface_quadrature_degree
ReferenceFiniteElements.surface_quadrature_points_and_weights
```

## Topology Interface
The below methods must be implemented to define the element 
toplogy.
```@docs
boundary_element
boundary_normals
dimension
edge_vertices
face_vertices
num_boundaries
num_edges
num_faces
num_vertices_per_cell
vertex_coordinates
```

# Element subtypes
Below are additional abstract types subtyped off of ```AbstractElementType``` for element
topologies of different dimensions and further subtyped based on common element topologies.

## 0-Dimensional Element types
Different types of 0-D elements (e.g. vertices of various implementations) can be implemented by subtyping off of the abstract type
```@docs
ReferenceFiniteElements.AbstractVertex
```

## 1-Dimensional Element types
Different types of 1-D elements (e.g. edges, sides, lines, etc. of various implementations) can be implemented by subtyping off of the abstract type
```@docs
ReferenceFiniteElements.AbstractEdge
```

## 2-Dimensional Element types
Different types of 2-D elements (e.g. faces, triangles, quads, polygons, etc. of various implementations) can be implemented by subtyping off of the abstract type
```@docs
ReferenceFiniteElements.AbstractFace
```
There are also various subtypes of this type including
```@docs
ReferenceFiniteElements.AbstractQuad
ReferenceFiniteElements.AbstractTri
```

## 3-Dimensional Element types
Different types of 3-D elements (e.g. volumes, hexes, tets, etc. of various implementations) can be implemented by subtyping off of the abstract type
```@docs
ReferenceFiniteElements.AbstractVolume
```
There are also various subtypes of this type including
```@docs
ReferenceFiniteElements.AbstractHex
ReferenceFiniteElements.AbstractPyramid
ReferenceFiniteElements.AbstractTet
ReferenceFiniteElements.AbstractWedge
```
