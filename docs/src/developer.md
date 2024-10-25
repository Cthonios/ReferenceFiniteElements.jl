```@meta
CurrentModule = ReferenceFiniteElements
```

# How to add a new element type
To add a new finite element there's a few things we need to do. First, we may need to define an abstract type if the element type is not a currently supported abstract element topology. The currently supported abstract element topologies are   
- ```AbstractEdge```
- ```AbstractHex```
- ```AbstractQuad```
- ```AbstractTet```
- ```AbstractTri```

We'll use a ```Tri3``` element as an example of how to add a new element type.

First, let's pretend ```AbstractTri``` doesn't exist yet. The first thing we need to do is define the following abstract type

```jldoctest developer_example
julia> using ReferenceFiniteElements

julia> abstract type AbstractTri{V, I, P, Q} <: ReferenceFiniteElements.AbstractFace{V, 3, I, P, Q} end

```
The four parameters ```V```, ```I```, ```P```, and ```Q``` are the number of vertices, interpolation type, the polynomial degree, and the number of quadrature degree for an element implementation. The dimension is fixed at ```2``` for the ```D``` parameter in the abstrac type. For a given element topology, such as a three noded triangle in this example, ```P``` will be fixed parameters while we allow for ```I``` and ```Q``` to be variable, i.e. have potentially different interpolation and quadrature rules on a 3 noded element. 

We now need to define some methods for basiscs about the element toplogy such
as number of vertices, edges, faces, etc.

e.g. 
```jldoctest developer_example
julia> num_edges(e::Type{<:AbstractTri}) = 3
num_edges (generic function with 1 method)


```

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


Now we define a struct for specific implementations. Here we'll define the Tri3 struct.
```jldoctest developer_example
julia> struct Tri3{I, Q} <: AbstractTri{3, I, 1, Q} end

julia> e = Tri3{Lagrange, 1}()
Tri3{Lagrange, 1}()

julia> num_edges(typeof(e))
3
```
Here we have fixed ```P``` to be ```1``` respectively while ```I``` and ```Q``` are left as parametric types.

We can also define some type stable constructors for different quadrature degrees. This look like the following

Now that we have a basic type to represent our element topology defined, there's a few methods we need to define for setup purposes. These methods can be split into methods need for abstract element topology types, such as ```AbstractTri```, and methods needed for specific element types such as ```Tri3```.


The methods needed for specific interpolation types include the following:
- ```shape_function_values```
- ```shape_funciton_gradients```
- ```shape_function_hessians```
