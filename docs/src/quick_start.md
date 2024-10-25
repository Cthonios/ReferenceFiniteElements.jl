# Quckstart
To setup a finite element you can do the following with a 8 node hexahedral element with a first order quadrature rule (8 integration points) used as an example below

```jldoctest quickstart
julia> using ReferenceFiniteElements

julia> re = ReferenceFE(Hex8{Lagrange, 2}())
ReferenceFE:
  Element Type                                  = Hex8{Lagrange, 2}
  Surface Element Type                          = Quad4{Lagrange, 2}
  Interpolation Type                            = Lagrange
  Polynomial Degree                             = 1
  Quadrature Degree                             = 2
  Array Backend Type                            = ReferenceFiniteElements.ArrayBackend{StaticArraysCore.SArray}()
  Edge Nodes Storage Type                       = Nothing
  Face Nodes Storage Type                       = Matrix{Int64}
  Interior Nodes Storage Type                   = Vector{Int64}
  Nodal Coordinates Storage Type                = Vector{StaticArraysCore.SVector{3, Float64}}
  CellInterpolants:
    Quadrature Point Storage Type               = Vector{StaticArraysCore.SVector{3, Float64}}
    Quadrature Weights Storage Type             = Vector{Float64}
    Shape Function Values Storage Type          = Vector{StaticArraysCore.SVector{8, Float64}}
    Shape Function Gradients Storage Type       = Vector{StaticArraysCore.SMatrix{8, 3, Float64, 24}}
    Shape Function Hessians Storage Type        = Vector{StaticArraysCore.SArray{Tuple{8, 3, 3}, Float64, 3, 72}}
  Surface Nodal Coordinates Storage Type        = Nothing
  CellInterpolants:
    Quadrature Point Storage Type               = Matrix{Vector{Float64}}
    Quadrature Weights Storage Type             = Matrix{Float64}
    Shape Function Values Storage Type          = Matrix{StaticArraysCore.SVector{8, Float64}}
    Shape Function Gradients Storage Type       = Matrix{StaticArraysCore.SMatrix{8, 3, Float64, 24}}
    Shape Function Hessians Storage Type        = Matrix{StaticArraysCore.SArray{Tuple{8, 3, 3}, Float64, 3, 72}}

```

To do something useful with our finite element we can look at the quadrature points. To get all the quadrature points, use the following method
```jldoctest quickstart
julia> quadrature_points(re)
8-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [-0.5773502691896258, -0.5773502691896258, -0.5773502691896258]
 [0.5773502691896258, -0.5773502691896258, -0.5773502691896258]
 [-0.5773502691896258, 0.5773502691896258, -0.5773502691896258]
 [0.5773502691896258, 0.5773502691896258, -0.5773502691896258]
 [-0.5773502691896258, -0.5773502691896258, 0.5773502691896258]
 [0.5773502691896258, -0.5773502691896258, 0.5773502691896258]
 [-0.5773502691896258, 0.5773502691896258, 0.5773502691896258]
 [0.5773502691896258, 0.5773502691896258, 0.5773502691896258]

```

To get a specific quadrature point, say the first point, we can use an analogous method as follows
```jldoctest quickstart
julia> quadrature_point(re, 1)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -0.5773502691896258
 -0.5773502691896258
 -0.5773502691896258

```

There are similar methods for the quadrature weights. See below.
```jldoctest quickstart
julia> quadrature_weights(re)
8-element Vector{Float64}:
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0

```

```jldoctest quickstart
julia> quadrature_weight(re, 1)
1.0

```

For the shape function values we can access these via the following methods
```jldoctest quickstart
julia> shape_function_values(re)
8-element Vector{StaticArraysCore.SVector{8, Float64}}:
 [0.4905626121623441, 0.13144585576580212, 0.035220810900864506, 0.13144585576580212, 0.13144585576580212, 0.035220810900864506, 0.009437387837655926, 0.035220810900864506]
 [0.13144585576580212, 0.4905626121623441, 0.13144585576580212, 0.035220810900864506, 0.035220810900864506, 0.13144585576580212, 0.035220810900864506, 0.009437387837655926]
 [0.13144585576580212, 0.035220810900864506, 0.13144585576580212, 0.4905626121623441, 0.035220810900864506, 0.009437387837655926, 0.035220810900864506, 0.13144585576580212]
 [0.035220810900864506, 0.13144585576580212, 0.4905626121623441, 0.13144585576580212, 0.009437387837655926, 0.035220810900864506, 0.13144585576580212, 0.035220810900864506]
 [0.13144585576580212, 0.035220810900864506, 0.009437387837655926, 0.035220810900864506, 0.4905626121623441, 0.13144585576580212, 0.035220810900864506, 0.13144585576580212]
 [0.035220810900864506, 0.13144585576580212, 0.035220810900864506, 0.009437387837655926, 0.13144585576580212, 0.4905626121623441, 0.13144585576580212, 0.035220810900864506]
 [0.035220810900864506, 0.009437387837655926, 0.035220810900864506, 0.13144585576580212, 0.13144585576580212, 0.035220810900864506, 0.13144585576580212, 0.4905626121623441]
 [0.009437387837655926, 0.035220810900864506, 0.13144585576580212, 0.035220810900864506, 0.035220810900864506, 0.13144585576580212, 0.4905626121623441, 0.13144585576580212]

```

```jldoctest quickstart
julia> shape_function_value(re, 1)
8-element StaticArraysCore.SVector{8, Float64} with indices SOneTo(8):
 0.4905626121623441
 0.13144585576580212
 0.035220810900864506
 0.13144585576580212
 0.13144585576580212
 0.035220810900864506
 0.009437387837655926
 0.035220810900864506

```

For shape function gradients, one can use the methods ```shape_function_gradients``` 
or ```shape_function_gradient``` respectively. This example is ommitted for brevity.

For shape function hessians, one can use the methods ```shape_function_hessians``` 
or ```shape_function_hessian``` respectively. This example is ommitted for brevity.
