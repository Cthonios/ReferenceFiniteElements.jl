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

function surface_quadrature_points_and_weights(e::E, t1::Type{A}, t2::Type{T}) where {
  A <: Union{SVector, MVector}, T <: Number, E <: AbstractTet
}
  D = num_dimensions(e)

  ξs_2d, ws_2d = quadrature_points_and_weights(Tri3(degree(e)), t1, t2)

  ξs = Matrix{A{D, T}}(undef, 4, length(ξs_2d))
  ws = Matrix{T}(undef, 4, length(ξs_2d))
  ns = Matrix{A{D, T}}(undef, 4, length(ξs_2d))

  # side 1
  for (q, ξ) in enumerate(ξs_2d)
    ξs[1, q] = A{D, T}(0.0, ξ[1], ξ[2])
    ws[1, q] = ws_2d[q]
    ns[1, q] = A{D, T}(-1., 0., 0.)
  end

  # side 2
  for (q, ξ) in enumerate(ξs_2d)
    ξs[2, q] = A{D, T}(ξ[1], ξ[2], 1. - ξ[1] - ξ[2])
    ws[2, q] = ws_2d[q]
    ns[2, q] = A{D, T}(1. / sqrt(3.), 1. / sqrt(3.), 1. / sqrt(3.))
  end

  # side 3
  for (q, ξ) in enumerate(ξs_2d)
    ξs[3, q] = A{D, T}(ξ[1], ξ[2], 0.0)
    ws[3, q] = ws_2d[q]
    ns[3, q] = A{D, T}(0., 0., -1.)
  end

  # side 4
  for (q, ξ) in enumerate(ξs_2d)
    ξs[4, q] = A{D, T}(ξ[1], 0.0, ξ[2])
    ws[4, q] = ws_2d[q]
    ns[4, q] = A{D, T}(0., -1., 0.)
  end

  return ξs, ws, ns
end