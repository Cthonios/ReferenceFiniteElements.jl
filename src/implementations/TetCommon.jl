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
struct Tet10 <: AbstractTet{10, 3}
  degree::Int64
end

types_to_generate = (
  (:(SVector{3, Float32}), :Float32),
  (:(SVector{3, Float64}), :Float64),
  (:(MVector{3, Float32}), :Float32),
  (:(MVector{3, Float64}), :Float64),
  # (:(Vector{Float64}), :Float64)
)

for type_pair in types_to_generate

  """
  """
  @eval function quadrature_points_and_weights(
    e::E, 
    ::Type{$(type_pair[1])},
    ::Type{$(type_pair[2])}
  ) where E <: AbstractTet

    if degree(e) == 1
      ξs    = Vector{$(type_pair[1])}(undef, 1)
      ξs[1] = $(type_pair[1])(1. / 4., 1. / 4., 1. / 4.)
      ws    = $(type_pair[2])[1. / 6.]
    elseif degree(e) == 2
      ξs    = Vector{$(type_pair[1])}(undef, 5)
      ξs[1] = $(type_pair[1])(1. / 4., 1. / 4., 1. / 4.)
      ξs[2] = $(type_pair[1])(1. / 6., 1. / 6., 1. / 6.)
      ξs[3] = $(type_pair[1])(1. / 6., 1. / 6., 1. / 2.)
      ξs[4] = $(type_pair[1])(1. / 6., 1. / 2., 1. / 6.)
      ξs[5] = $(type_pair[1])(1. / 2., 1. / 6., 1. / 6.)

      #
      ws = $(type_pair[2])[
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

end
