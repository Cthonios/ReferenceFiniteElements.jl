# ReferenceFiniteElements 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cthonios.github.io/ReferenceFiniteElements.jl/) 
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cthonios.github.io/ReferenceFiniteElements.jl/dev/) 
[![Build Status](https://github.com/Cthonios/ReferenceFiniteElements.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Cthonios/ReferenceFiniteElements.jl/actions/workflows/CI.yml?query=branch%3Amain) 
[![Coverage](https://app.codecov.io/gh/Cthonios/ReferenceFiniteElements.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Cthonios/ReferenceFiniteElements.jl)

ReferenceFiniteElements.jl is meant to serve as an educational and lightweight package to define finite elements in the reference (quadrature space) configuration. The goal is to provide a very simple to utilize interface so researchers can play with new finite element formulations without having to learn an entire large code or package base. 

```
julia> re = ReferenceFE(Quad4(2))
Element type             = ReferenceFE{Int64, 4, 2, Float64, 8, 16, StructArrays.StructVector{ReferenceFiniteElements.Interpolants{4, 2, Float64, 8, 16}, NamedTuple{(:ξ, :w, :N, :∇N_ξ, :∇∇N_ξ), Tuple{Vector{StaticArraysCore.SVector{2, Float64}}, Vector{Float64}, Vector{StaticArraysCore.SVector{4, Float64}}, Vector{StaticArraysCore.SMatrix{4, 2, Float64, 8}}, Vector{StaticArraysCore.SArray{Tuple{4, 2, 2}, Float64, 3, 16}}}}, Int64}, Matrix{Float64}, Matrix{Int64}, Vector{Int64}}

Nodal coordinates        = 
2×4 Matrix{Float64}:
 -1.0   1.0  1.0  -1.0
 -1.0  -1.0  1.0   1.0
Face nodes               = 
2×4 Matrix{Int64}:
 1  2  3  4
 2  3  4  1

Int64[]

Shape function values    = 
4-element StaticArraysCore.SVector{4, Float64} with indices SOneTo(4):
 0.6220084679281462
 0.16666666666666663
 0.044658198738520435
 0.16666666666666663
4-element StaticArraysCore.SVector{4, Float64} with indices SOneTo(4):
 0.16666666666666663
 0.6220084679281462
 0.16666666666666663
 0.044658198738520435
4-element StaticArraysCore.SVector{4, Float64} with indices SOneTo(4):
 0.16666666666666663
 0.044658198738520435
 0.16666666666666663
 0.6220084679281462
4-element StaticArraysCore.SVector{4, Float64} with indices SOneTo(4):
 0.044658198738520435
 0.16666666666666663
 0.6220084679281462
 0.16666666666666663

Shape function gradients = 
4×2 StaticArraysCore.SMatrix{4, 2, Float64, 8} with indices SOneTo(4)×SOneTo(2):
 -0.394338  -0.394338
  0.394338  -0.105662
  0.105662   0.105662
 -0.105662   0.394338
4×2 StaticArraysCore.SMatrix{4, 2, Float64, 8} with indices SOneTo(4)×SOneTo(2):
 -0.394338  -0.105662
  0.394338  -0.394338
  0.105662   0.394338
 -0.105662   0.105662
4×2 StaticArraysCore.SMatrix{4, 2, Float64, 8} with indices SOneTo(4)×SOneTo(2):
 -0.105662  -0.394338
  0.105662  -0.105662
  0.394338   0.105662
 -0.394338   0.394338
4×2 StaticArraysCore.SMatrix{4, 2, Float64, 8} with indices SOneTo(4)×SOneTo(2):
 -0.105662  -0.105662
  0.105662  -0.394338
  0.394338   0.394338
 -0.394338   0.105662

Shape function hessians  = 
4×2×2 StaticArraysCore.SArray{Tuple{4, 2, 2}, Float64, 3, 16} with indices SOneTo(4)×SOneTo(2)×SOneTo(2):
[:, :, 1] =
 0.0   0.25
 0.0  -0.25
 0.0   0.25
 0.0  -0.25

[:, :, 2] =
  0.25  0.0
 -0.25  0.0
  0.25  0.0
 -0.25  0.0
4×2×2 StaticArraysCore.SArray{Tuple{4, 2, 2}, Float64, 3, 16} with indices SOneTo(4)×SOneTo(2)×SOneTo(2):
[:, :, 1] =
 0.0   0.25
 0.0  -0.25
 0.0   0.25
 0.0  -0.25

[:, :, 2] =
  0.25  0.0
 -0.25  0.0
  0.25  0.0
 -0.25  0.0
4×2×2 StaticArraysCore.SArray{Tuple{4, 2, 2}, Float64, 3, 16} with indices SOneTo(4)×SOneTo(2)×SOneTo(2):
[:, :, 1] =
 0.0   0.25
 0.0  -0.25
 0.0   0.25
 0.0  -0.25

[:, :, 2] =
  0.25  0.0
 -0.25  0.0
  0.25  0.0
 -0.25  0.0
4×2×2 StaticArraysCore.SArray{Tuple{4, 2, 2}, Float64, 3, 16} with indices SOneTo(4)×SOneTo(2)×SOneTo(2):
[:, :, 1] =
 0.0   0.25
 0.0  -0.25
 0.0   0.25
 0.0  -0.25

[:, :, 2] =
  0.25  0.0
 -0.25  0.0
  0.25  0.0
 -0.25  0.0
```
