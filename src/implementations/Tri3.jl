"""
"""
function ReferenceFEStencil(::Tri3, degree::I, Itype::Type = Integer, Rtype::Type = Float64) where I <: Integer
  points::Matrix{Rtype} = [
    0.0 1.0 0.0;
    0.0 0.0 1.0
  ]
  vertex_points::Vector{Itype} = [1, 2, 3]
  face_points::Matrix{Itype} = [
    1 2 3
    2 3 1
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return ReferenceFEStencil{Itype, Rtype, Tri3}(degree, points, vertex_points, face_points, interior_nodes)
end

function shape_function_values_int(::Tri3, ξ)
  N = @SVector [
    1. - ξ[1] - ξ[2],
    ξ[1],
    ξ[2]
  ]
  # N = MVector{3, Float64}(1. - ξ[1] - ξ[2], ξ[1], ξ[2])
end

function shape_function_gradients_int(::Tri3, ξ)
  ∇N_ξ = @SMatrix [
    -1. -1.;
     1.  0.;
     0.  1.
  ]
end
