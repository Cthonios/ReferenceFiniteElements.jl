"""
"""
function ReferenceFEStencil(::Quad4, degree::I, Itype::Type = Integer, Rtype::Type = Float64) where I <: Integer
  points = [
    -1.0  1.0 1.0 -1.0;
    -1.0 -1.0 1.0  1.0
  ]
  vertex_points = [1, 2, 3, 4]
  face_points = [
    1 2 3 4
    2 3 4 1
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return ReferenceFEStencil{Itype, Rtype, Quad4}(degree, points, vertex_points, face_points, interior_nodes)
end

function shape_function_values_int(::Quad4, ξ)
  N = @SVector [
    0.25 * (1.0 - ξ[1]) * (1.0 - ξ[2]),
    0.25 * (1.0 + ξ[1]) * (1.0 - ξ[2]),
    0.25 * (1.0 + ξ[1]) * (1.0 + ξ[2]),
    0.25 * (1.0 - ξ[1]) * (1.0 + ξ[2])
  ]
end

function shape_function_gradients_int(::Quad4, ξ)
  ∇N_ξ = @SMatrix [
    -0.25 * (1.0 - ξ[2]) -0.25 * (1.0 - ξ[1]);
     0.25 * (1.0 - ξ[2]) -0.25 * (1.0 + ξ[1]);
     0.25 * (1.0 + ξ[2])  0.25 * (1.0 + ξ[1]);
    -0.25 * (1.0 + ξ[2])  0.25 * (1.0 - ξ[1])
  ]
end
