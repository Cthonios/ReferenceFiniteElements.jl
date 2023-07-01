"""
"""
struct Tri6 <: AbstractTri
end

"""
"""
function ElementStencil(e::Tri6, degree::I, Itype::Type = Integer, Rtype::Type = Float64) where I <: Integer
  points = [
    0.0 1.0 0.0 0.5 0.5 0.0;
    0.0 0.0 1.0 0.0 0.5 0.5
  ]
  vertex_points = [1, 2, 3]
  face_points = [
    1 2 3
    4 5 6
    2 3 1
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return ElementStencil{Itype, Rtype, Tri6}(e, degree, points, vertex_points, face_points, interior_nodes)
end

function shape_function_values(::Tri6, ξ)
  λ = 1. - ξ[1] - ξ[2]
  N = @SVector [
    λ * (2. * λ - 1.),
    ξ[1] * (2. * ξ[1] - 1.),
    ξ[2] * (2. * ξ[2] - 1.),
    4. * ξ[1] * λ,
    4. * ξ[1] * ξ[2],
    4. * ξ[2] * λ
  ]
end

function shape_function_gradients(::Tri6, ξ)
  λ = 1. - ξ[1] - ξ[2]
  ∇N_ξ = @SMatrix [
    -1. * (2. * λ - 1.) - 2. * λ       -1. * (2. * λ - 1.) - 2. * λ      ;
     (2. * ξ[1] - 1.) + 2. * ξ[1]       0.0                              ;
     0.0                                (2. * ξ[2] - 1.) + 2. * ξ[2]     ;
     4. * λ - 4. * ξ[1]                -4. * ξ[1]                        ;
     4. * ξ[2]                          4. * ξ[1]                        ;
    -4. * ξ[2]                          4. * λ - 4. * ξ[2]               ;  
  ]
end 

export Tri6
