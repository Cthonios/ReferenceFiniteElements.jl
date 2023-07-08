"""
"""
function element_stencil(::Tet10, degree::I, ::Type{Itype}, ::Type{Ftype}) where {I <: Integer, Itype <: Integer, Ftype <: AbstractFloat}
  # first one stupidly causes an allocation
  # nodal_coordinates = Ftype[
  #   0.0 1.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0
  #   0.0 0.0 1.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5
  #   0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.5 0.5 0.5
  # ]
  nodal_coordinates = Matrix{Ftype}(undef, 3, 10)
  nodal_coordinates[1, 1] = 0.0
  nodal_coordinates[2, 1] = 0.0
  nodal_coordinates[3, 1] = 0.0
  #
  nodal_coordinates[1, 2] = 1.0
  nodal_coordinates[2, 2] = 0.0
  nodal_coordinates[3, 2] = 0.0
  #
  nodal_coordinates[1, 3] = 0.0
  nodal_coordinates[2, 3] = 1.0
  nodal_coordinates[3, 3] = 0.0
  #
  nodal_coordinates[1, 4] = 0.0
  nodal_coordinates[2, 4] = 0.0
  nodal_coordinates[3, 4] = 1.0
  #
  nodal_coordinates[1, 5] = 0.5
  nodal_coordinates[2, 5] = 0.0
  nodal_coordinates[3, 5] = 0.0
  #
  nodal_coordinates[1, 6] = 0.5
  nodal_coordinates[2, 6] = 0.5
  nodal_coordinates[3, 6] = 0.0
  #
  nodal_coordinates[1, 7] = 0.0
  nodal_coordinates[2, 7] = 0.5
  nodal_coordinates[3, 7] = 0.0
  #
  nodal_coordinates[1, 8] = 0.0
  nodal_coordinates[2, 8] = 0.0
  nodal_coordinates[3, 8] = 0.5
  #
  nodal_coordinates[1, 9] = 0.5
  nodal_coordinates[2, 9] = 0.0
  nodal_coordinates[3, 9] = 0.5
  #
  nodal_coordinates[1, 10] = 0.0
  nodal_coordinates[2, 10] = 0.5
  nodal_coordinates[3, 10] = 0.5
  #

  # first one stupidly causes an allocation
  # face_nodes = Itype[
  #   1 2  1  1
  #   5 6  8  7
  #   2 3  4  3
  #   9 10 10 6
  #   4 4  3  2
  #   8 8  7  5
  # ]
  face_nodes = Matrix{Itype}(undef, 6, 4)
  face_nodes[1, 1] = 1
  face_nodes[2, 1] = 5
  face_nodes[3, 1] = 2
  face_nodes[4, 1] = 9
  face_nodes[5, 1] = 4
  face_nodes[6, 1] = 8
  #
  face_nodes[1, 2] = 2
  face_nodes[2, 2] = 6
  face_nodes[3, 2] = 3
  face_nodes[4, 2] = 10
  face_nodes[5, 2] = 4
  face_nodes[6, 2] = 8
  #
  face_nodes[1, 3] = 1
  face_nodes[2, 3] = 8
  face_nodes[3, 3] = 4
  face_nodes[4, 3] = 10
  face_nodes[5, 3] = 3
  face_nodes[6, 3] = 7
  #
  face_nodes[1, 4] = 1
  face_nodes[2, 4] = 7
  face_nodes[3, 4] = 3
  face_nodes[4, 4] = 6
  face_nodes[5, 4] = 2
  face_nodes[6, 4] = 5
  #
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, face_nodes, interior_nodes
end

function shape_function_values(::Tet10, ξ::T) where T <: AbstractArray
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

function shape_function_gradients(::Tet10, ξ::T) where T <: AbstractArray
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
