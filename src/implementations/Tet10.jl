"""
"""
struct Tet10 <: AbstractTet
end

"""
"""
function ElementStencil(e::Tet10, degree::I, Itype::Type = Integer, Rtype::Type = Float64) where I <: Integer
  points = [
    0.0 1.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0
    0.0 0.0 1.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5
    0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.5 0.5 0.5
  ]
  vertex_points = [1, 2, 3, 4]
  face_points = [
    1 2  1  1
    5 6  8  7
    2 3  4  3
    9 10 10 6
    4 4  3  2
    8 8  7  5
  ]
  interior_nodes = []
  return ElementStencil{Itype, Rtype}(e, degree, points, vertex_points, face_points, interior_nodes)
end

function shape_function_values(::Tet10, ξ)
  λ = 1.0 - ξ[1] - ξ[2] - ξ[3]
  N = @SVector [
    λ * (2.0 * λ - 1.0),
    ξ[1] * (2.0 * ξ[1] - 1.0),
    ξ[2] * (2.0 * ξ[2] - 1.0),
    ξ[3] * (2.0 * ξ[3] - 1.0),
    4.0 * ξ[1] * λ,
    4.0 * ξ[1] * ξ[2],
    4.0 * ξ[2] * λ,
    4.0 * ξ[3] * λ,
    4.0 * ξ[1] * ξ[3],
    4.0 * ξ[2] * ξ[3]
  ]
end

function shape_function_gradients(::Tet10, ξ)
  λ = 1.0 - ξ[1] - ξ[2] - ξ[3]
  ∇N_ξ = @SMatrix [
    -(2.0 * λ - 1.0) - 2.0 * λ -(2.0 * λ - 1.0) - 2.0 * λ  -(2.0 * λ - 1.0) - 2.0 * λ;
    (2.0 * ξ[1] - 1.0) + 2.0 * ξ[1] 0.0 0.0;
    0.0 (2.0 * ξ[2] - 1.0) + 2.0 * ξ[2] 0.0;
    0.0 0.0 (2.0 * ξ[3] - 1.0) + 2.0 * ξ[3];
    4.0 * λ - 4.0 * ξ[1] -4.0 * ξ[1] -4.0 * ξ[1];
    4.0 * ξ[2] 4.0 * ξ[1] 0.0;
    -4.0 * ξ[2] 4.0 * λ - 4.0 * ξ[2] 0.0;
    -4.0 * ξ[3] -4.0 * ξ[3] 4.0 * λ - 4.0 * ξ[3];
    4.0 * ξ[3] 0.0 4.0 * ξ[1];
    0.0 4.0 * ξ[3] 4.0 * ξ[2]
  ]
end

export Tet10