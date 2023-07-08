"""
"""
function element_stencil(::Quad9, degree::I, ::Type{Itype}, ::Type{Ftype}) where {I <: Integer, Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates = Ftype[
    -1.0  1.0 1.0 -1.0  0.0 1.0 0.0 -1.0 0.0
    -1.0 -1.0 1.0  1.0 -1.0 0.0 1.0  0.0 0.0
  ]
  # @time face_nodes = Itype[
  #   1 2 3 4
  #   5 6 7 8
  #   2 3 4 1
  # ]
  # allocation if I don't do the below rather than above
  face_nodes = Matrix{Itype}(undef, 3, 4)
  face_nodes[1, 1] = 1
  face_nodes[1, 2] = 2
  face_nodes[1, 3] = 3
  face_nodes[1, 4] = 4
  #
  face_nodes[2, 1] = 5
  face_nodes[2, 2] = 6
  face_nodes[2, 3] = 7
  face_nodes[2, 4] = 8
  #
  face_nodes[3, 1] = 2
  face_nodes[3, 2] = 3
  face_nodes[3, 3] = 4
  face_nodes[3, 4] = 1

  interior_nodes = Itype[9]
  return nodal_coordinates, face_nodes, interior_nodes
end

function shape_function_values(::Quad9, ξ::T) where T <: AbstractArray
  N = @SVector [
    0.25 * (ξ[1]^2 - ξ[1]) * (ξ[2]^2 - ξ[2]),
    0.25 * (ξ[1]^2 + ξ[1]) * (ξ[2]^2 - ξ[2]),
    0.25 * (ξ[1]^2 + ξ[1]) * (ξ[2]^2 + ξ[2]),
    0.25 * (ξ[1]^2 - ξ[1]) * (ξ[2]^2 + ξ[2]),
    0.50 * (ξ[2]^2 - ξ[2]) * (1.0 - ξ[1]^2),
    0.50 * (ξ[1]^2 + ξ[1]) * (1.0 - ξ[2]^2), 
    0.50 * (ξ[2]^2 + ξ[2]) * (1.0 - ξ[1]^2), 
    0.50 * (ξ[1]^2 - ξ[1]) * (1.0 - ξ[2]^2), 
    (1.0 - ξ[1]^2) * (1.0 - ξ[2]^2)
  ]
end

function shape_function_gradients(::Quad9, ξ::T) where T <: AbstractArray
  ∇N_ξ = @SMatrix [
    0.25 * (ξ[2]^2 - ξ[2]) * (2. * ξ[1] - 1.)  0.25 * (ξ[1]^2 - ξ[1]) * (2. * ξ[2] - 1.);
    0.25 * (ξ[2]^2 - ξ[2]) * (2. * ξ[1] + 1.)  0.25 * (ξ[1]^2 + ξ[1]) * (2. * ξ[2] - 1.);
    0.25 * (ξ[2]^2 + ξ[2]) * (2. * ξ[1] + 1.)  0.25 * (ξ[1]^2 + ξ[1]) * (2. * ξ[2] + 1.);
    0.25 * (ξ[2]^2 + ξ[2]) * (2. * ξ[1] - 1.)  0.25 * (ξ[1]^2 - ξ[1]) * (2. * ξ[2] + 1.);
    0.50 * (ξ[2]^2 - ξ[2]) * (-2.0 * ξ[1])     0.50 * (2.0 * ξ[2] - 1.0) * (1.0 - ξ[1]^2);
    0.50 * (2.0 * ξ[1] + 1.0) * (1.0 - ξ[2]^2) 0.50 * (ξ[1]^2 + ξ[1]) * (-2.0 * ξ[2]);
    0.50 * (ξ[2]^2 + ξ[2]) * (-2.0 * ξ[1])     0.50 * (2.0 * ξ[2] + 1.0) * (1.0 - ξ[1]^2);
    0.50 * (2.0 * ξ[1] - 1.0) * (1.0 - ξ[2]^2) 0.50 * (ξ[1]^2 - ξ[1]) * (-2.0 * ξ[2]);
    (-2.0 * ξ[1]) * (1.0 - ξ[2]^2)             (1.0 - ξ[1]^2) * (-2.0 * ξ[2])
  ]
end
