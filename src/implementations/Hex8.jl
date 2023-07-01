"""
"""
struct Hex8 <: AbstractHex
end

"""
"""
function ElementStencil(e::Hex8, degree::I, Itype::Type = Integer, Rtype::Type = Float64) where I <: Integer
  points = [
    -1.0  1.0  1.0 -1.0 -1.0  1.0 1.0 -1.0
    -1.0 -1.0  1.0  1.0 -1.0 -1.0 1.0  1.0
    -1.0 -1.0 -1.0  -1.0 1.0  1.0 1.0  1.0
  ]
  vertex_points = [1, 2, 3, 4, 5, 6, 7, 8]
  face_points = [
    1 2 3 1 1 5
    2 3 4 5 4 6
    6 7 8 8 3 7
    5 6 7 4 2 8
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return ElementStencil{Itype, Rtype, Hex8}(e, degree, points, vertex_points, face_points, interior_nodes)
end

function shape_function_values(::Hex8, ξ)
  N = @SVector [
    0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]) * (1.0 - ξ[3]),
    0.125 * (1.0 + ξ[1]) * (1.0 - ξ[2]) * (1.0 - ξ[3]),
    0.125 * (1.0 + ξ[1]) * (1.0 + ξ[2]) * (1.0 - ξ[3]),
    0.125 * (1.0 - ξ[1]) * (1.0 + ξ[2]) * (1.0 - ξ[3]),
    0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]) * (1.0 + ξ[3]),
    0.125 * (1.0 + ξ[1]) * (1.0 - ξ[2]) * (1.0 + ξ[3]),
    0.125 * (1.0 + ξ[1]) * (1.0 + ξ[2]) * (1.0 + ξ[3]),
    0.125 * (1.0 - ξ[1]) * (1.0 + ξ[2]) * (1.0 + ξ[3]),
  ]
end

function shape_function_gradients(::Hex8, ξ)
  ∇N_ξ = @SMatrix [
    -0.125 * (1.0 - ξ[2]) * (1.0 - ξ[3]) -0.125 * (1.0 - ξ[1]) * (1.0 - ξ[3]) -0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]);
     0.125 * (1.0 - ξ[2]) * (1.0 - ξ[3]) -0.125 * (1.0 + ξ[1]) * (1.0 - ξ[3]) -0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]);
     0.125 * (1.0 + ξ[2]) * (1.0 - ξ[3])  0.125 * (1.0 + ξ[1]) * (1.0 - ξ[3]) -0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]);
    -0.125 * (1.0 + ξ[2]) * (1.0 - ξ[3])  0.125 * (1.0 - ξ[1]) * (1.0 - ξ[3]) -0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]);
    -0.125 * (1.0 - ξ[2]) * (1.0 + ξ[3]) -0.125 * (1.0 - ξ[1]) * (1.0 + ξ[3])  0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]);
     0.125 * (1.0 - ξ[2]) * (1.0 + ξ[3]) -0.125 * (1.0 + ξ[1]) * (1.0 + ξ[3])  0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]);
     0.125 * (1.0 + ξ[2]) * (1.0 + ξ[3])  0.125 * (1.0 + ξ[1]) * (1.0 + ξ[3])  0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]);
    -0.125 * (1.0 + ξ[2]) * (1.0 + ξ[3])  0.125 * (1.0 - ξ[1]) * (1.0 + ξ[3])  0.125 * (1.0 - ξ[1]) * (1.0 - ξ[2]);
  ]
end

export Hex8
