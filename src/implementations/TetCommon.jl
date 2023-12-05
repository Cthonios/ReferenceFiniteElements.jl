"""
"""
abstract type AbstractTet{N, D, Q} <: ReferenceFEType{N, D, Q} end

"""
"""
struct Tet4{Q} <: AbstractTet{4, 3, Q}
  degree::Int64
end

"""
"""
struct Tet10{Q} <: AbstractTet{10, 3, Q}
  degree::Int64
end

Tet4(::Val{1})  = Tet4{1}(1)
Tet4(::Val{2})  = Tet4{5}(2)
Tet10(::Val{1}) = Tet10{1}(1)
Tet10(::Val{2}) = Tet10{5}(2)

# using Tet4(1) as template
for type_pair in types_to_generate_quadrature(Tet4(Val(1)))

  """
  """
  @eval function quadrature_points_and_weights(
    e::E, 
    ::Type{$(type_pair[2])},
    ::Type{$(type_pair[1])}
  ) where E <: AbstractTet

    if degree(e) == 1
      ξs    = Vector{$(type_pair[2])}(undef, 1)
      ξs[1] = $(type_pair[2])(1. / 4., 1. / 4., 1. / 4.)
      ws    = $(type_pair[1])[1. / 6.]
    elseif degree(e) == 2
      ξs    = Vector{$(type_pair[2])}(undef, 5)
      ξs[1] = $(type_pair[2])(1. / 4., 1. / 4., 1. / 4.)
      ξs[2] = $(type_pair[2])(1. / 6., 1. / 6., 1. / 6.)
      ξs[3] = $(type_pair[2])(1. / 6., 1. / 6., 1. / 2.)
      ξs[4] = $(type_pair[2])(1. / 6., 1. / 2., 1. / 6.)
      ξs[5] = $(type_pair[2])(1. / 2., 1. / 6., 1. / 6.)

      #
      ws = $(type_pair[1])[
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
