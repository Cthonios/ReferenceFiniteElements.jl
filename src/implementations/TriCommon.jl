"""
"""
abstract type AbstractTri{N, D, Q} <: ReferenceFEType{N, D, Q} end

"""
"""
struct Tri3{Q} <: AbstractTri{3, 2, Q}
  degree::Int64
end

Tri3(::Val{1}) = Tri3{1}(1)
Tri3(::Val{2}) = Tri3{3}(2)
Tri3(::Val{3}) = Tri3{6}(3)
Tri3(::Val{4}) = Tri3{6}(4)
Tri3(::Val{5}) = Tri3{7}(5)
Tri3(::Val{6}) = Tri3{12}(6)
Tri3(q::Int)   = Tri3(Val{q}())

"""
"""
struct Tri6{Q} <: AbstractTri{6, 2, Q}
  degree::Int64
end

Tri6(::Val{1}) = Tri6{1}(1)
Tri6(::Val{2}) = Tri6{3}(2)
Tri6(::Val{3}) = Tri6{6}(3)
Tri6(::Val{4}) = Tri6{6}(4)
Tri6(::Val{5}) = Tri6{7}(5)
Tri6(::Val{6}) = Tri6{12}(6)
Tri6(q::Int)   = Tri6(Val{q}())

"""
"""
struct SimplexTri{N, Q} <: AbstractTri{N, 2, Q}
  n_nodes::Int64
  degree::Int64
end

# SimplexTri(::Val{1}, ::Val{1}) = SimplexTri{3, 1}(3, 1)
# SimplexTri(::Val{1}, ::Val{2}) = SimplexTri{3, 3}(3, 2)
# SimplexTri(::Val{1}, ::Val{3}) = SimplexTri{3, 6}(3, 3)
# SimplexTri(::Val{1}, ::Val{4}) = SimplexTri{3, 6}(3, 4)
# SimplexTri(::Val{1}, ::Val{5}) = SimplexTri{3, 7}(3, 5)
# SimplexTri(::Val{1}, ::Val{6}) = SimplexTri{3, 12}(3, 6)

# SimplexTri(::Val{2}, ::Val{1}) = SimplexTri{6, 1}(6, 1)
# SimplexTri(::Val{2}, ::Val{2}) = SimplexTri{6, 3}(6, 2)
# SimplexTri(::Val{2}, ::Val{3}) = SimplexTri{6, 6}(6, 3)
# SimplexTri(::Val{2}, ::Val{4}) = SimplexTri{6, 6}(6, 4)
# SimplexTri(::Val{2}, ::Val{5}) = SimplexTri{6, 7}(6, 5)
# SimplexTri(::Val{2}, ::Val{6}) = SimplexTri{6, 12}(6, 6)

# SimplexTri(n::Int, q::Int)     = SimplexTri(Val(n), Val(q))

SimplexTri(::Val{1}) = SimplexTri{3, 1}(3, 1)
SimplexTri(::Val{2}) = SimplexTri{6, 3}(6, 2)
SimplexTri(::Val{3}) = SimplexTri{10, 6}(10, 3)
SimplexTri(q::Int)   = SimplexTri(Val(q))

function quadrature_points_and_weights(e::E, ::Type{A}, ::Type{T}) where {
  A <: Union{SVector, MVector}, T <: Number, E <: AbstractTri
}

  D = num_dimensions(e)

  if degree(e) == 1
    ξs    = Vector{A{D, T}}(undef, 1)
    ξs[1] = A{D, T}(1. / 3., 1. / 3.)
    ws    = T[0.5]
  elseif degree(e) == 2
    ξs    = Vector{A{D, T}}(undef, 3)
    ξs[1] = A{D, T}(2. / 3., 1. / 6.)
    ξs[2] = A{D, T}(1. / 6., 2. / 3.)
    ξs[3] = A{D, T}(1. / 6., 1. / 6.)
    ws    = T[1. / 6., 1. / 6., 1. / 6.]
  elseif degree(e) <= 4
    ξs    = Vector{A{D, T}}(undef, 6)
    ξs[1] = A{D, T}(1.081030181680700E-01, 4.459484909159650E-01)
    ξs[2] = A{D, T}(4.459484909159650E-01, 1.081030181680700E-01)
    ξs[3] = A{D, T}(4.459484909159650E-01, 4.459484909159650E-01)
    ξs[4] = A{D, T}(8.168475729804590E-01, 9.157621350977100E-02)
    ξs[5] = A{D, T}(9.157621350977100E-02, 8.168475729804590E-01)
    ξs[6] = A{D, T}(9.157621350977100E-02, 9.157621350977100E-02)

    ws    = T[
      1.116907948390055E-01,
      1.116907948390055E-01,
      1.116907948390055E-01,
      5.497587182766100E-02,
      5.497587182766100E-02,
      5.497587182766100E-02
    ]
  elseif degree(e) <= 5
    ξs    = Vector{A{D, T}}(undef, 7)
    ξs[1] = A{D, T}(3.33333333333333E-01, 3.33333333333333E-01)
    ξs[2] = A{D, T}(5.97158717897700E-02, 4.70142064105115E-01)
    ξs[3] = A{D, T}(4.70142064105115E-01, 5.97158717897700E-02)
    ξs[4] = A{D, T}(4.70142064105115E-01, 4.70142064105115E-01)
    ξs[5] = A{D, T}(7.97426985353087E-01, 1.01286507323456E-01)
    ξs[6] = A{D, T}(1.01286507323456E-01, 7.97426985353087E-01)
    ξs[7] = A{D, T}(1.01286507323456E-01, 1.01286507323456E-01)

    ws    = T[
      1.12500000000000E-01,
      6.61970763942530E-02,
      6.61970763942530E-02,
      6.61970763942530E-02,
      6.29695902724135E-02,
      6.29695902724135E-02,
      6.29695902724135E-02
    ]
  elseif degree(e) <= 6
    ξs     = Vector{A{D, T}}(undef, 12)
    ξs[1]  = A{D, T}(5.01426509658179E-01, 2.49286745170910E-01)
    ξs[2]  = A{D, T}(2.49286745170910E-01, 5.01426509658179E-01)
    ξs[3]  = A{D, T}(2.49286745170910E-01, 2.49286745170910E-01)
    ξs[4]  = A{D, T}(8.73821971016996E-01, 6.30890144915020E-02)
    ξs[5]  = A{D, T}(6.30890144915020E-02, 8.73821971016996E-01)
    ξs[6]  = A{D, T}(6.30890144915020E-02, 6.30890144915020E-02)
    ξs[7]  = A{D, T}(5.31450498448170E-02, 3.10352451033784E-01)
    ξs[8]  = A{D, T}(6.36502499121399E-01, 5.31450498448170E-02)
    ξs[9]  = A{D, T}(3.10352451033784E-01, 6.36502499121399E-01)
    ξs[10] = A{D, T}(5.31450498448170E-02, 6.36502499121399E-01)
    ξs[11] = A{D, T}(6.36502499121399E-01, 3.10352451033784E-01)
    ξs[12] = A{D, T}(3.10352451033784E-01, 5.31450498448170E-02)

    ws     = T[
      5.83931378631895E-02,
      5.83931378631895E-02,
      5.83931378631895E-02,
      2.54224531851035E-02,
      2.54224531851035E-02,
      2.54224531851035E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02
    ]
  end
  return ξs, ws
end

function surface_quadrature_points_and_weights(e::E, ::Type{A}, ::Type{T}) where {
  A <: Union{SVector, MVector}, T <: Number, E <: AbstractTri
}
  D = num_dimensions(e)
  ξs_1d, ws_1d = gausslegendre(degree(e))
  ξs_1d .= (ξs_1d .+ 1.) ./ 2.
  ws_1d .= 0.5 * ws_1d

  ξs = Matrix{A{D, T}}(undef, 3, length(ξs_1d))
  ws = Matrix{T}(undef, 3, length(ξs_1d))
  ns = Matrix{A{D, T}}(undef, 3, length(ξs_1d))

  # side 1
  for (q, ξ) in enumerate(ξs_1d)
    ξs[1, q] = A{D, T}(ξ, 0.0)
    ws[1, q] = ws_1d[q]
    ns[1, q] = A{D, T}(ξ, -1.0)
  end

  # side 2
  for (q, ξ) in enumerate(ξs_1d)
    ξs[2, q] = A{D, T}(ξ, 1 - ξ)
    ws[2, q] = ws_1d[q]
    ns[2, q] = A{D, T}(1. / sqrt(2.), 1. / sqrt(2.))
  end

  # side 3
  for (q, ξ) in enumerate(ξs_1d)
    ξs[3, q] = A{D, T}(0.0, ξ)
    ws[3, q] = ws_1d[q]
    ns[3, q] = A{D, T}(-1.0, ξ)
  end

  return ξs, ws, ns
end