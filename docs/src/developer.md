```@meta
CurrentModule = ReferenceFiniteElements
```

# How to add a new element type
To add a new finite element there's a few things we need to do. First, we may need to define an abstract type if the element type is not a currently supported abstract element topology. The currently supported abstract element topologies are   
- ```AbstractEdge```
- ```AbstractFace```
- ```AbstractQuad```
- ```MyAbstractTri```
- ```AbstractVolume```
- ```AbstractHex```
- ```AbstractTet```


We'll use a ```MyTri3``` element as an example of how to add a new element type.

First, let's pretend ```MyAbstractTri``` doesn't exist yet. The first thing we need to do is define the following abstract type

```jldoctest developer_example
julia> using ReferenceFiniteElements

julia> abstract type MyAbstractTri{PT, PD} <: ReferenceFiniteElements.AbstractFace{PT, PD} end

```
The two parameters ```PT``` and ```PD``` are the polynomial type and polynomial degree.

There are four loose interfaces associated with each element type
- The topology interace: this defines methods that describes basics of an element such as vertices, edges, etc.
- The dof interface: this defines methods describes where/how dofs are used on an element for a given interpolation type
- The quadrature interface: this defines methods that implement quadrature stencils
- The shape function interface: this implement interpolants on the element

## Topology interface
```julia
function boundary_element end
function boundary_normals end # always returns 3 x num boundaries
function dimension end
function edge_vertices end
function face_vertices end
function num_boundaries end
function num_edges end
function num_faces end
function num_vertices_per_cell end
function vertex_coordinates end # always returns 3 x num_vertices_per_cell
```
For a face, the ```dimension```, ```face_vertices```, ```num_boundaries```, and ```num_faces``` methods are already implemented since they are constant for any face type.

Now we define an abstract type for a triangle
```jldoctest developer_example
julia> abstract type MyAbstractTri{PT, PD} <: ReferenceFiniteElements.AbstractFace{PT, PD} end

```
We can finish the topology interface on this abstract type with the following method implementations
```jldoctest developer_example
julia> boundary_element(::MyAbstractTri{PT, PD}, ::Int) where {PT, PD} = Edge{PT, PD}(; shifted = true)
boundary_element (generic function with 1 method)

julia> boundary_normals(::MyAbstractTri) = [0. 1. / sqrt(2.) -1.; -1. 1. / sqrt(2.) 0.; 0. 0. 0.]
boundary_normals (generic function with 1 method)

julia> edge_vertices(::MyAbstractTri) = [1 2 3; 2 3 1]
edge_vertices (generic function with 1 method)

julia> num_edges(::MyAbstractTri) = 3
num_edges (generic function with 1 method)

julia> num_vertices_per_cell(::MyAbstractTri) = 3
num_vertices_per_cell (generic function with 1 method)

julia> vertex_coordinates(::MyAbstractTri) = [0. 1. 0.; 0. 0. 1.; 0. 0. 0.]
vertex_coordinates (generic function with 1 method)
```

## Dof interface
Now we define a struct for specific implementation. In this case the implementation is general. In the future we could have specific triangle implementations that break this pattern such as ```Tri7``` or complex composite type elements. Here we'll define the ```MyTri``` struct.
```jldoctest developer_example
julia> struct MyTri{PD, PT} <: MyAbstractTri{PD, PT} end
```


## Quadrature interface
Below is a simple implementation of the quadrature stencil for first and second order.
```jldoctest developer_example
function cell_quadrature_points_and_weights(::MyAbstractTri, q_rule::GaussLobattoLegendre)
  if cell_quadrature_degree(q_rule) == 1
    ξs = Matrix{Float64}(undef, 2, 1)
    ξs[:, 1] = [1. / 3., 1. / 3.]
    ws = [0.5]
  elseif cell_quadrature_degree(q_rule) == 2
    ξs = Matrix{Float64}(undef, 2, 3)
    ξs[:, 1] = [2. / 3., 1. / 6.]
    ξs[:, 2] = [1. / 6., 2. / 3.]
    ξs[:, 3] = [1. / 6., 1. / 6.]
    ws = [1. / 6., 1. / 6., 1. / 6.]
  else
    @assert false
  end
end

# output
cell_quadrature_points_and_weights (generic function with 1 method)
```

## Shape function interface
Now that we have a basic type to represent our element topology defined, we now need to 
implement our interpolation scheme. The methods needed for specific interpolation types include the following:
- ```shape_function_values```
- ```shape_funciton_gradients```
- ```shape_function_hessians```

Below is an example for the ```MyTri3``` element above
```jldoctest developer_example
function shape_function_value(e::MyTri{Lagrange, PD}, X, ξ) where PD
  Ns = [
    1. - ξ[1] - ξ[2],
    ξ[1],
    ξ[2]
  ]
  return Ns
end

# output
shape_function_value (generic function with 1 method)

```

```jldoctest developer_example
function shape_function_gradient(e::MyTri{Lagrange, PD}, X, ξ) where PD
  Ns = [
    -1. -1.;
     1.  0.; 
     0.  1.
  ]
  return Ns
end

# output
shape_function_gradient (generic function with 1 method)

```

```jldoctest developer_example
function shape_function_hessian(e::MyTri{Lagrange, PD}, X, ξ) where PD
  Ns = zeros(3, 2, 2)
  return Ns
end

# output
shape_function_hessian (generic function with 1 method)

```
