"""
"""
function ReferenceFEStencil(::Tet4, degree::I, Itype::Type = Integer, Rtype::Type = Float64) where I <: Integer
  points = [
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 1.0
  ]
  vertex_points = [1, 2, 3, 4]
  face_points = [
    1 2 1 1
    2 3 4 3
    4 4 3 2
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return ReferenceFEStencil{Itype, Rtype, Tet4}(degree, points, vertex_points, face_points, interior_nodes)
end

function shape_function_values_int(::Tet4, ξ)
  N = @SVector [
    1. - ξ[1] - ξ[2] - ξ[3],
    ξ[1],
    ξ[2],
    ξ[3]
  ]
end

function shape_function_gradients_int(::Tet4, ξ)
  ∇N_ξ = @SMatrix [
    -1.0 1.0 0.0 0.0
    -1.0 0.0 1.0 0.0
    -1.0 0.0 0.0 1.0
  ]
end
