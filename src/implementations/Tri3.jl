"""
"""
struct Tri3 <: AbstractTri
end

"""
"""
function ElementStencil(e::Tri3, degree::I, Itype::Type = Integer, Rtype::Type = Float64) where I <: Integer
  points = [
    0.0 1.0 0.0;
    0.0 0.0 1.0
  ]
  vertex_points = [1, 2, 3]
  face_points = [
    1 2 3
    2 3 1
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return ElementStencil{Itype, Rtype, Tri3}(e, degree, points, vertex_points, face_points, interior_nodes)
end

function shape_function_values(::Tri3, ξ)
  N = @SVector [
    1. - ξ[1] - ξ[2],
    ξ[1],
    ξ[2]
  ]
end

function shape_function_gradients(::Tri3, ξ)
  ∇N_ξ = @SMatrix [
    -1. -1.;
     1.  0.;
     0.  1.
  ]
end

export Tri3
