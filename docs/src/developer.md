```@meta
CurrentModule = ReferenceFiniteElements
```

# How to add a new element type
To add a new finite element there's a few things we need to do. First, we may need to define an abstract type if the element type is not a currently supported abstract element topology. The currently supported abstract element topologies are   
- `AbstractEdge`
- `AbstractHex`
- `AbstractQuad`
- `AbstractTet`
- `AbstractTri`

We'll use a `Tri3` element as an example of how to add a new element type.

First, let's pretend `AbstractTri` doesn't exist yet. The first thing we need to do is define the following abstract type

```julia
abstract type AbstractTri{N, D, Q} <: RefenceFEType{N, D, Q} end
```
The three parameters `N`, `D`, and `Q` are the number of nodes, the number of dimensions, and the number of quadrature points for an element implementation. For a given element topology, such as a three noded triangle in this example, `N` and `D` will be fixed parameters while we allow for `Q` to be variable, i.e. have potentially different quadrature rules on the element. 

Now we define a struct for specific implementations. Here we'll define the Tri3 struct.
```julia
struct Tri3{Q} <: AbstractTri{3, 2, Q}
  degree::Int
end
```
Here we have fixed `N` and `D` to be `3` and `2` respectively while `Q` is left as a parametric type. The struct has one field, `degree`, which denotes the degree of the quadrature rule. Other element rules could have more parameters depending upon how interpolation and integration are performed. 

We can also define some type stable constructors for different quadrature degrees. This look like the following

```julia
Tri3(::Val{1}) = Tri3{1}(1)
Tri3(::Val{2}) = Tri3{3}(2)
```
for quadrature degrees of `1` and `2` respectively. The numbers `1` and `3` correspond to the number of quadrature points in each rule.

Now that we have a basic type to represent our element topology defined, there's a few methods we need to define for setup purposes. These methods are

- `quadrature_points_and_weights`
- `element_stencil`
- `shape_function_values`
- `shape_funciton_gradients`
- `shape_function_hessians`

Below is the implementaiton of `element_stencial` for `AbstractTri` elements

```julia
function quadrature_points_and_weights(e::E, ::Type{A}, ::Type{T}) where {
  A <: Union{SVector, MVector}, T <: Number, E <: AbstractTri
}

  D = num_dimensions(e)

  if degree(e) == 1
    ξs    = Vector{A{D, T}}(undef, 1)
    ξs[1] = A{D, T}(1. / 3., 1. / 3.)
    ws    = T[0.5]
  elseif degree(e) == 2
    ξs    = Vector{A{D, T}}(undef, 3)
    ξs[1] = A{D, T}(2. / 3., 1. / 6.)
    ξs[2] = A{D, T}(1. / 6., 2. / 3.)
    ξs[3] = A{D, T}(1. / 6., 1. / 6.)
    ws    = T[1. / 6., 1. / 6., 1. / 6.]
  end
  return ξs, ws
end
```
Two different quadature rules, `degree = 1` and `degree = 2` are supported. The `degree = 1` case we could refer to as "fully integrated."

The implementation of `element_stencil` for `Tri3` elements is shown below

```julia
function element_stencil(::Tri3, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates = Ftype[
    0.0 1.0 0.0;
    0.0 0.0 1.0
  ]
  face_nodes = Itype[
    1 2 3
    2 3 1
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, face_nodes, interior_nodes
end
```

The implementation of `shape_function_values` for `Tri3` elements is shown below

```julia
function shape_function_values(e::Tri3, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SVector, MVector}, A2 <: AbstractArray{<:Number, 1}
}
  N = A1{num_nodes(e), eltype(ξ)}(
    1. - ξ[1] - ξ[2],
    ξ[1],
    ξ[2]
  )
end
```
The implementation of `shape_function_gradients` for `Tri3` elements is shown below

```julia
function shape_function_gradients(e::Tri3, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SMatrix, MMatrix}, A2 <: AbstractArray{<:Number, 1}
}
  ∇N_ξ = A1{num_nodes(e), num_dimensions(e), eltype(ξ), num_nodes(e) * num_dimensions(e)}(
    -1., 1., 0.,
    -1., 0., 1.
  )
end
```
The implementation of `shape_function_hessians` for `Tri3` elements is shown below
```julia
function shape_function_hessians(e::Tri3, ::Type{A1}, ξ::A2) where {
  A1 <: Union{SArray, MArray}, A2 <: AbstractArray{<:Number, 1}
}
  N, D = num_nodes(e), num_dimensions(e)
  ∇∇N_ξ = zeros(A1{Tuple{N, D, D}, eltype(ξ), 3, N * D * D})
end
```

As a convention, all of the interpolation methods for every element type are defined with `StaticArrays` as the default array type. Internally, when other array types are requested, these arrays are converted to the requested type. This allows for allocations on setup only for those array types that require it. 

Another thing to note is for all existing elements defined in the package, the `quadrature_points_and_weights` method is defined on abstract types. This could change in the future, but is the current convention. 
