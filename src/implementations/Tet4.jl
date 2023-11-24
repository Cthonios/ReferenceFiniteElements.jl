"""
"""
function element_stencil(::Tet4, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  # first one stupidly causes an allocation
  # nodal_coordinates = Ftype[
  #   0.0 1.0 0.0 0.0
  #   0.0 0.0 1.0 0.0
  #   0.0 0.0 0.0 1.0
  # ]
  nodal_coordinates = Matrix{Ftype}(undef, 3, 4)
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
  # first one stupidly causes an allocation
  # face_nodes = Itype[
  #   1 2 1 1
  #   2 3 4 3
  #   4 4 3 2
  # ]
  face_nodes = Matrix{Itype}(undef, 3, 4)
  face_nodes[1, 1] = 1
  face_nodes[2, 1] = 2
  face_nodes[3, 1] = 4
  #
  face_nodes[1, 2] = 2
  face_nodes[2, 2] = 3
  face_nodes[3, 2] = 4
  #
  face_nodes[1, 3] = 1
  face_nodes[2, 3] = 4
  face_nodes[3, 3] = 3
  #
  face_nodes[1, 4] = 1
  face_nodes[2, 4] = 3
  face_nodes[3, 4] = 2
  #
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, face_nodes, interior_nodes
end

# using Tet4(1) as a template
for type in types_to_generate_interpolants(Tet4(1))
  """
  """
  @eval function shape_function_values(e::Tet4, ::Type{$(type[1])}, ξ::A) where A <: AbstractArray{<:Number, 1}
    N = $(type[1]){num_nodes(e), eltype(ξ)}(
      1. - ξ[1] - ξ[2] - ξ[3],
      ξ[1],
      ξ[2],
      ξ[3]
    )
  end
end

"""
"""
function shape_function_gradients(::Tet4, ξ::SVector{3, <:Real})
  ∇N_ξ = (@SMatrix [
    -1. -1. -1.;
     1.  0.  0.;
     0.  1.  0.;
     0.  0.  1.;
  ]) |> SMatrix{4, 3, eltype(ξ), 12}
end

"""
"""
function shape_function_hessians(::Tet4, ξ::SVector{3, <:Real})
  ∇∇N_ξ = zeros(SArray{Tuple{4, 3, 3}, eltype(ξ), 3, 36})
end