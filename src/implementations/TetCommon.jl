"""
"""
const Tet4 = ReferenceFEType{4, 3}
const Tet10 = ReferenceFEType{10, 3}
const TetUnion = Union{Tet4, Tet10}

"""
"""
function Quadrature(::E, degree::I, Rtype::Type = Float64) where {E <: TetUnion, I <: Integer}
  if degree == 1
    ξ = Matrix{Rtype}(undef, 3, 1)
    ξ[1, 1] = 1. / 4.
    ξ[2, 1] = 1. / 4.
    ξ[3, 1] = 1. / 4.
    w = [1. / 6.]
  elseif degree == 2
    ξ = [
      1. / 4. 1. / 6. 1. / 6. 1. / 6. 1. / 2. 
      1. / 4. 1. / 6. 1. / 6. 1. / 2. 1. / 6.
      1. / 4. 1. / 6. 1. / 2. 1. / 6. 1. / 6.
    ]
    w = [
      -2. / 15.
      3. / 40.
      3. / 40.
      3. / 40.
      3. / 40.
    ]
  else
    throw(ErrorException("Unsupported quadrature egree"))
  end
  return Quadrature{Rtype}(ξ, w)
end
