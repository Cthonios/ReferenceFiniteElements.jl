# ReferenceFiniteElements [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cmhamel.github.io/ReferenceFiniteElements.jl/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cmhamel.github.io/ReferenceFiniteElements.jl/dev/) [![Build Status](https://github.com/cmhamel/ReferenceFiniteElements.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cmhamel/ReferenceFiniteElements.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/cmhamel/ReferenceFiniteElements.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Cthonios/ReferenceFiniteElements.jl)

ReferenceFiniteElements.jl is meant to serve as an educational and lightweight package to define finite elements in the reference (quadrature space) configuration. The goal is to provide a very simple to utilize interface so researchers can play with new finite element formulations without having to learn an entire large code or package base. 

```
julia> re = ReferenceFE(Quad4(), 1)
ReferenceFE{4, 2, Integer, Float64, Quad4}(Quadrature:
  Quadrature points = [0.0, 0.0]
                      
  Quadrature weights = 4.0
                       
, ParentElement:
  Coordinates: [-1.0, -1.0]
               [1.0, -1.0]
               [1.0, 1.0]
               [-1.0, 1.0]
               
  Vertex nodes: Integer[1, 2, 3, 4]

  Face nodes: Integer[1, 2]
              Integer[2, 3]
              Integer[3, 4]
              Integer[4, 1]
              
  Interior nodes: Integer[]
, ShapeFunctionPair{4, 2, Float64}[ShapeFunctionPair{4, 2, Float64}([0.25, 0.25, 0.25, 0.25], [-0.25 -0.25; 0.25 -0.25; 0.25 0.25; -0.25 0.25])])
```
