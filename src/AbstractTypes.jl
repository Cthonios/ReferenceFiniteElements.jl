# D - dimension
# P - polynomial
# Q - quadrature order
# T - type of interpolation
abstract type AbstractElementType{I, D, P, Q} end
dimension(::AbstractElementType{I, D, P, Q}) where {I, D, P, Q} = D
interpolation_type(::AbstractElementType{I, D, P, Q}) where {I, D, P, Q} = I
polynomial_degree(::AbstractElementType{I, D, P, Q}) where {I, D, P, Q} = P
quadrature_degree(::AbstractElementType{I, D, P, Q}) where {I, D, P, Q} = Q

function num_surfaces(e::AbstractElementType)
  if dimension(e) == 0
    return 0
  elseif dimension(e) == 1
    return 2
  elseif dimension(e) == 2
    return num_edges(e)
  elseif dimension(e) == 3
    return num_faces(e)
  end
end

abstract type AbstractInterpolationType end

abstract type AbstractInterpolantsContainer{I} end
