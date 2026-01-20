# Quckstart
To setup a finite element you can do the following with a 8 node hexahedral element with a first order quadrature rule (8 integration points) used as an example below

```jldoctest quickstart
julia> using ReferenceFiniteElements

julia> re = ReferenceFE(Hex{Lagrange, 1}(), GaussLobattoLegendre(2))
ReferenceFE
  Element type      = Hex{Lagrange, 1}()
  Polynomial type   = Lagrange
  Polynomial degree = 1



```

To do something useful with our finite element we can look at the quadrature points. To get a specific quadrature point, say the first point, we can use an analogous method as follows
```jldoctest quickstart
julia> cell_quadrature_point(re, 1)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -0.5773502691896258
 -0.5773502691896258
 -0.5773502691896258

```

There are similar methods for the quadrature weights. See below.
```jldoctest quickstart
julia> cell_quadrature_weight(re, 1)
1.0

```

For the shape function values we can access these via the following methods
```jldoctest quickstart
julia> cell_shape_function_value(re, 1)
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
