# type used to exploit multiple dispatch
"""
Type to define new element shapes
"""
abstract type ReferenceFEType{N, D, P, Q} end

"""
Returns the quadrature degree of a ReferenceFEType
"""
degree(e::ReferenceFEType) = e.degree

"""
Returns the number of nodes for a ReferenceFEType
"""
num_nodes(::ReferenceFEType{N, D, P, Q}) where {N, D, P, Q} = N

"""
Returns the number of dimensions for a ReferenceFEType
"""
num_dimensions(::ReferenceFEType{N, D, P, Q}) where {N, D, P, Q} = D

"""
Polynomial degree, i.e. p-refinement
"""
polynomial_degree(::ReferenceFEType{N, D, P, Q}) where {N, D, P, Q} = P

"""
Returns the number of quadrature porints for a ReferenceFEType
"""
num_q_points(::ReferenceFEType{N, D, P, Q}) where {N, D, P, Q} = Q

# functions to be defined
function quadrature_points_and_weights end
function shape_function_values end
function shape_function_gradients end
function shape_function_hessians end

function quadrature_points_and_weights(e::ReferenceFEType, ::Type{Vector}, T::Type{<:Number})
  ξs, ws = quadrature_points_and_weights(e, SVector, T)
  new_ξs = map(ξ -> Vector(ξ), ξs)
  return new_ξs, ws
end

function shape_function_values(e::ReferenceFEType, ::Type{Vector}, ξ::A) where A <: AbstractArray{<:Number, 1}
  Ns = shape_function_values(e, SVector, ξ)
  new_Ns = map(N -> Vector(N), Ns)
  return new_Ns
end

function shape_function_gradients(e::ReferenceFEType, ::Type{Matrix}, ξ::A) where A <: AbstractArray{<:Number, 1}
  ∇N_ξs     = shape_function_gradients(e, SMatrix, ξ)
  new_∇N_ξs = map(∇N_ξ -> Matrix(∇N_ξ), ∇N_ξs)
  return new_∇N_ξs
end

function shape_function_hessians(e::ReferenceFEType, ::Type{Array}, ξ::A) where A <: AbstractArray{<:Number, 1}
  ∇∇N_ξs     = shape_function_hessians(e, SArray, ξ)
  new_∇∇N_ξs = map(∇∇N_ξ -> Array(∇∇N_ξ), ∇∇N_ξs)
  return new_∇∇N_ξs
end