"""
"""
abstract type AbstractTet{N, D} <: ReferenceFEType{N, D} end

"""
"""
struct Tet4 <: AbstractTet{4, 3}
  degree::Int64
end

"""
"""
function quadrature_points_and_weights(e::E, ::Type{Ftype}) where {E <: AbstractTet, Ftype <: AbstractFloat}
  if degree(e) == 1
    # ξs = Matrix{Ftype}(undef, 3, 1)
    # ξs[1, 1] = 1. / 4.
    # ξs[2, 1] = 1. / 4.
    # ξs[3, 1] = 1. / 4.
    ξs = Vector{SVector{3, Ftype}}(undef, 1)
    ξs[1] = SVector{3, Ftype}(1. / 4., 1. / 4., 1. / 4.)
    ws = Ftype[1. / 6.]
  elseif degree(e) == 2
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
