module ReferenceFiniteElementsForwardDiffExt

export AutoDiff

using ForwardDiff
using ReferenceFiniteElements
using StaticArrays

struct AutoDiff
end

ReferenceFiniteElements.shape_function_gradients(
  ::Type{AutoDiff}, ::R, ξ::S
) where {R <: ReferenceFiniteElements.ReferenceFEType, S <: SVector} = 
ForwardDiff.jacobian(x -> shape_function_values(e, x), ξ)

end # module