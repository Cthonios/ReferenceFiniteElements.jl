"""
"""
function element_stencil(::Quad4, degree::I, ::Type{Itype}, ::Type{Ftype}) where {I <: Integer, Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates = Ftype[
    -1.0  1.0 1.0 -1.0;
    -1.0 -1.0 1.0  1.0
  ]
  face_nodes = Itype[
    1 2 3 4
    2 3 4 1
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, face_nodes, interior_nodes
end

function shape_function_values(::Quad4, ξ::SVector{2, Ftype}) where Ftype <: AbstractFloat
  N = SVector{4, Ftype}(
    0.25 * (1.0 - ξ[1]) * (1.0 - ξ[2]),
    0.25 * (1.0 + ξ[1]) * (1.0 - ξ[2]),
    0.25 * (1.0 + ξ[1]) * (1.0 + ξ[2]),
    0.25 * (1.0 - ξ[1]) * (1.0 + ξ[2])
  )
end

function shape_function_gradients(::Quad4, ξ::SVector{2, Ftype}) where Ftype <: AbstractFloat
  # ∇N_ξ = @SMatrix [
  #   -0.25 * (1.0 - ξ[2]) -0.25 * (1.0 - ξ[1]);
  #    0.25 * (1.0 - ξ[2]) -0.25 * (1.0 + ξ[1]);
  #    0.25 * (1.0 + ξ[2])  0.25 * (1.0 + ξ[1]);
  #   -0.25 * (1.0 + ξ[2])  0.25 * (1.0 - ξ[1])
  # ]
  ∇N_ξ = SMatrix{4, 2, Ftype, 8}(
    -0.25 * (1.0 - ξ[2]), 0.25 * (1.0 - ξ[2]), 0.25 * (1.0 + ξ[2]), -0.25 * (1.0 + ξ[2]),
    #
    -0.25 * (1.0 - ξ[1]), -0.25 * (1.0 + ξ[1]), 0.25 * (1.0 + ξ[1]), 0.25 * (1.0 - ξ[1])
  )
end
