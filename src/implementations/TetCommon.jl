"""
"""
const Tet4 = ReferenceFEType{4, 3}
const Tet10 = ReferenceFEType{10, 3}
const TetUnion = Union{Tet4, Tet10}

"""
"""
function quadrature_points_and_weights(::E, degree::I, ::Type{Ftype} = Float64) where {E <: TetUnion, I <: Integer, Ftype <: AbstractFloat}
  if degree == 1
    ξs = Matrix{Ftype}(undef, 3, 1)
    ξs[1, 1] = 1. / 4.
    ξs[2, 1] = 1. / 4.
    ξs[3, 1] = 1. / 4.
    ws = Ftype[1. / 6.]
  elseif degree == 2
    # @time ξs = [
    #   1. / 4. 1. / 6. 1. / 6. 1. / 6. 1. / 2. 
    #   1. / 4. 1. / 6. 1. / 6. 1. / 2. 1. / 6.
    #   1. / 4. 1. / 6. 1. / 2. 1. / 6. 1. / 6.
    # ]
    # allocation if I don't do below
    ξs = Matrix{Ftype}(undef, 3, 5)
    ξs[1, 1] = 1. / 4.
    ξs[2, 1] = 1. / 4.
    ξs[3, 1] = 1. / 4.
    #
    ξs[1, 2] = 1. / 6.
    ξs[2, 2] = 1. / 6.
    ξs[3, 2] = 1. / 6.
    #
    ξs[1, 3] = 1. / 6.
    ξs[2, 3] = 1. / 6.
    ξs[3, 3] = 1. / 2.
    #
    ξs[1, 4] = 1. / 6.
    ξs[2, 4] = 1. / 2.
    ξs[3, 4] = 1. / 6.
    #
    ξs[1, 5] = 1. / 2.
    ξs[2, 5] = 1. / 6.
    ξs[3, 5] = 1. / 6.
    #
    ws = Ftype[
      -2. / 15.
       3. / 40.
       3. / 40.
       3. / 40.
       3. / 40.
    ]
  else
    throw(ErrorException("Unsupported quadrature egree"))
  end
  return ξs, ws
end
