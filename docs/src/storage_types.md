<!-- ```@meta
CurrentModule = ReferenceFiniteElements
```

# Storage types
There are different storage types available for the main type `ReferenceFE` beyond the default settings of `StructArrays` composed of `StaticArrays`. To use these default types we can use the following optional keyword arguments in the constructors as follows

```jldoctest
using ReferenceFiniteElements
using StaticArrays

re = ReferenceFE(
  Hex8(2); 
  int_type=Int64, float_type=Float64, 
  array_type=MArray, storage_type=Array
)

# output

ReferenceFE
  Element type                = Hex8{8}
  Dimension                   = 3
  Number of nodes             = 8
  Number of quadrature points = 8
  Integer type                = Int64
  Float type                  = Float64
  Nodal coordinates type      = Matrix{Float64}
  Edge nodes type             = Matrix{Int64}
  Face nodes type             = Matrix{Int64}
  Interior nodes type         = Vector{Int64}
  Interpolants type           = Interpolants{MArray, MArray, MArray}

```
The example above now uses mutable `MArray`s rather than the default of `SArray`. The storage type is also now a `Vector` of `MArray`s rather than a `StructArray` of `MArray`s. The currently supported internal array types are `Array`, `MArray`, and `SArray`. The currently supported storage types are `Array` and `StructArray`. -->
