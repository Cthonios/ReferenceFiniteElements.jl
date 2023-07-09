"""
"""
const Tet4 = ReferenceFEType{4, 3}
const Tet10 = ReferenceFEType{10, 3}
const TetUnion = Union{Tet4, Tet10}

"""
"""
function quadrature_points_and_weights(::E, degree::I, ::Type{Ftype}) where {E <: TetUnion, I <: Integer, Ftype <: AbstractFloat}
  if degree == 1
    # ξs = Matrix{Ftype}(undef, 3, 1)
    # ξs[1, 1] = 1. / 4.
    # ξs[2, 1] = 1. / 4.
    # ξs[3, 1] = 1. / 4.
    ξs = Vector{SVector{3, Ftype}}(undef, 1)
    ξs[1] = SVector{3, Ftype}(1. / 4., 1. / 4., 1. / 4.)
    ws = Ftype[1. / 6.]
  elseif degree == 2
    # @time ξs = [
    #   1. / 4. 1. / 6. 1. / 6. 1. / 6. 1. / 2. 
    #   1. / 4. 1. / 6. 1. / 6. 1. / 2. 1. / 6.
    #   1. / 4. 1. / 6. 1. / 2. 1. / 6. 1. / 6.
    # ]
    ξs = Vector{SVector{3, Ftype}}(undef, 5)
    ξs[1] = SVector{3, Ftype}(1. / 4., 1. / 4., 1. / 4.)
    ξs[2] = SVector{3, Ftype}(1. / 6., 1. / 6., 1. / 6.)
    ξs[3] = SVector{3, Ftype}(1. / 6., 1. / 6., 1. / 2.)
    ξs[4] = SVector{3, Ftype}(1. / 6., 1. / 2., 1. / 6.)
    ξs[5] = SVector{3, Ftype}(1. / 2., 1. / 6., 1. / 6.)

    #
    ws = Ftype[
      -2. / 15.
       3. / 40.
       3. / 40.
       3. / 40.
       3. / 40.
    ]
  else
    throw(ErrorException("Unsupported quadrature degree"))
  end
  return ξs, ws
end
