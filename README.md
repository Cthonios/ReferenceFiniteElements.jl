# ReferenceFiniteElements 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cthonios.github.io/ReferenceFiniteElements.jl/) 
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cthonios.github.io/ReferenceFiniteElements.jl/dev/) 
[![Build Status](https://github.com/Cthonios/ReferenceFiniteElements.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Cthonios/ReferenceFiniteElements.jl/actions/workflows/CI.yml?query=branch%3Amain) 
[![Coverage](https://codecov.io/gh/Cthonios/ReferenceFiniteElements.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Cthonios/ReferenceFiniteElements.jl)

ReferenceFiniteElements.jl is meant to serve as an educational and lightweight package to define finite elements in the reference (quadrature space) configuration. The goal is to provide a very simple to utilize interface so researchers can play with new finite element formulations without having to learn an entire large code or package base. 

# Installation
To use ReferenceFiniteElements first install it via the package manager via the following command

```julia
pkg> add ReferenceFiniteElements
```

# Quckstart
To setup a finite element you can do the following with a 8 node hexahedral element with a second order quadrature rule (8 integration points) used as an example below

```julia
using ReferenceFiniteElements
re = ReferenceFE(Hex8{Lagrange, 2}())

# output

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

```julia
quadrature_points(re)
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

```julia
ξ = quadrature_point(re, 1)

# output

3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -0.5773502691896258
 -0.5773502691896258
 -0.5773502691896258

```

There are similar methods for the quadrature weights. See below.

```julia
ws = quadrature_weights(re)

# output

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

```julia
w = quadrature_weights(re, 1)

# output

1.0

```

For the shape function values we can access these via the following methods

```julia
Ns = shape_function_values(re)

# output

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

```julia
N = shape_function_values(re, 1)

# output

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
