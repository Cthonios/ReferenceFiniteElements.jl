```@meta
CurrentModule = ReferenceFiniteElements
```

# How to add a new element type
To add a new finite element there's a few things we need to do. First, we may need to define an abstract type if the element type is not a currently supported abstract element topology. The currently supported abstract element topologies are   
- ```AbstractEdge```
- ```AbstractFace```
- ```AbstractQuad```
- ```AbstractTri```
- ```AbstractVolume```
- ```AbstractHex```
- ```AbstractTet```


We'll use a ```MyTri3``` element as an example of how to add a new element type.

First, let's pretend ```AbstractTri``` doesn't exist yet. The first thing we need to do is define the following abstract type

```jldoctest developer_example
julia> using ReferenceFiniteElements

julia> abstract type AbstractTri{V, I, P, Q} <: ReferenceFiniteElements.AbstractFace{V, 3, I, P, Q} end

```
The four parameters ```V```, ```I```, ```P```, and ```Q``` are the number of vertices, interpolation type, the polynomial degree, and the number of quadrature degree for an element implementation. For a given element topology, such as a three noded triangle in this example, ```V``` will be fixed parameters while we allow for ```I```, ```P```, and ```Q``` to be variable, i.e. have potentially different interpolation and quadrature rules on a triangular element. The number 3 above corresponds to the number of edges in the face topology. 

We now need to define some methods for basics about the element toplogy such
as number of vertices, edges, faces, etc.

TODO finish this part of the documention!
TODO change ```nodal_coordinates``` and ```surface_nodal_coordinates``` to be
```vertex_coordinates``` and ```surface_vertex_coordinates``` respectively.
TODO also ```surface_nodal_coordinates``` is likely redundant.

The complete list of methods that need to be defined for a new abstract topology type include the following:
- ```element_edge_nodes```
- ```element_face_nodes```
- ```element_interior_nodes```
- ```nodal_coordinates```
- ```num_edges```
- ```num_faces```
- ```num_interior_vertices```
- ```num_quadrature_points```
- ```num_vertices```
- ```num_vertices_per_edge```
- ```num_vertices_per_face```
- ```quadrature_points_and_weights```
- ```surface_element```
- ```surface_element_type```
- ```surface_nodal_coordinates```
- ```surface_quadrature_points_and_weights```


Now we define a struct for specific implementations. Here we'll define the MyTri3 struct.
```jldoctest developer_example
julia> struct MyTri3{I, Q} <: AbstractTri{3, I, 1, Q} end

julia> ReferenceFiniteElements.surface_element(e::MyTri3{Lagrange, Q}) where Q = Edge2{Lagrange, Q}()

julia> e = MyTri3{Lagrange, 1}()
MyTri3{Lagrange, 1}()

julia> ReferenceFiniteElements.num_edges(e)
3
```
Here we have fixed ```P``` to be ```1``` respectively while ```I``` and ```Q``` are left as parametric types to correspond to a linear interpolation but with a generic interpolation and quadrature rule.

Now that we have a basic type to represent our element topology defined, we now need to 
implement our interpolation scheme. The methods needed for specific interpolation types include the following:
- ```shape_function_values```
- ```shape_funciton_gradients```
- ```shape_function_hessians```

Below is an example for the ```MyTri3``` element above
```jldoctest developer_example
function shape_function_value(e::MyTri3{Lagrange}, X, ξ, backend::ReferenceFiniteElements.ArrayBackend)
  Ns = convert_to_vector(e, backend,
    1. - ξ[1] - ξ[2],
    ξ[1],
    ξ[2]
  )
  return Ns
end

# output
shape_function_value (generic function with 1 method)

```

```jldoctest developer_example
function shape_function_gradient(e::MyTri3{Lagrange}, X, ξ, backend::ReferenceFiniteElements.ArrayBackend)
  Ns = convert_to_matrix(e, backend,
    -1., 
    1., 
    0.,
    #
    -1., 
    0., 
    1.
  )
  return Ns
end

# output
shape_function_gradient (generic function with 1 method)

```

```jldoctest developer_example
function shape_function_hessian(e::MyTri3{Lagrange}, X, ξ, backend::ReferenceFiniteElements.ArrayBackend)
  Ns = convert_to_3d_array(e, backend,
    0., 0.,
    0., 0.,
    0., 0.,
    #
    0., 0.,
    0., 0.,
    0., 0.
  )
  return Ns
end

# output
shape_function_hessian (generic function with 1 method)

```
