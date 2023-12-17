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

Tet4(q::Int64)  = Tet4(Val{q}())
Tet10(q::Int64) = Tet10(Val{q}())

function quadrature_points_and_weights(e::E, ::Type{A}, ::Type{T}) where {
  E <: AbstractTet, A <: Union{SVector, MVector}, T <: Number
}

  D = num_dimensions(e)
  if degree(e) == 1
    ξs    = Vector{A{D, T}}(undef, 1)
    ξs[1] = A{D, T}(1. / 4., 1. / 4., 1. / 4.)
    ws    = T[1. / 6.]
  elseif degree(e) == 2
    ξs    = Vector{A{D, T}}(undef, 5)
    ξs[1] = A{D, T}(1. / 4., 1. / 4., 1. / 4.)
    ξs[2] = A{D, T}(1. / 6., 1. / 6., 1. / 6.)
    ξs[3] = A{D, T}(1. / 6., 1. / 6., 1. / 2.)
    ξs[4] = A{D, T}(1. / 6., 1. / 2., 1. / 6.)
    ξs[5] = A{D, T}(1. / 2., 1. / 6., 1. / 6.)

    #
    ws = T[
      -2. / 15.
       3. / 40.
       3. / 40.
       3. / 40.
       3. / 40.
    ]
  end
  return ξs, ws
end
