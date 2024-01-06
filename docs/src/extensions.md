# ReferenceFiniteElementsAdaptExt
For using `ReferenceFiniteElements` on GPUs. So far only tested with `CUDA.jl` as the GPU backend.

# ReferenceFiniteElementsExodusExt
Convenience methods for setting up `ReferenceFE`s based on an Exodus block. Really only useful if you use exodusII mesh types.

# ReferenceFiniteElementsRecipesBaseExt
Plotting tools for educational purposes. To use, try the following
```julia
using ReferenceFiniteElements, Plots, LaTeXStrings
re = ReferenceFE(Tri6(2))
plot(re) # plots nodes and edges of element along with quadrature points
plot(re, 1) # plots shape function values. The index is the shape function index
plot(re, 1, 1) # plots shape function gradients. The first index is for the index, second is for dimension
plot(re, 1, 1, 1) # plots shape function hessians. The first index  is for the shape function index and the second and third are for dimensions
```
