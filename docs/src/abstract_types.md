# Type hierarchy
The root of the type hierarchy tree in ```ReferenceFiniteElements``` is the following 
abstract type.

```@docs
ReferenceFiniteElements.AbstractElementType
```

There are some useful methods associated with this abstract type such as
```@docs
ReferenceFiniteElements.dimension
ReferenceFiniteElements.num_vertices
ReferenceFiniteElements.num_edges
ReferenceFiniteElements.num_faces
ReferenceFiniteElements.interpolation_type
ReferenceFiniteElements.polynomial_degree
ReferenceFiniteElements.quadrature_degree
```

That can fetch the parameteric types.

Other useful methods associated with this abstract type are the following which can be used to query information about the element topology.
```@docs
ReferenceFiniteElements.num_interior_vertices
ReferenceFiniteElements.num_quadrature_points
ReferenceFiniteElements.num_shape_functions
ReferenceFiniteElements.num_vertices_per_edge
ReferenceFiniteElements.num_vertices_per_face
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
ReferenceFiniteElements.AbstractTet
```

# Abstract Types for Storage Containers
```@docs
ReferenceFiniteElements.AbstractInterpolationType
ReferenceFiniteElements.AbstractInterpolantsContainer
```