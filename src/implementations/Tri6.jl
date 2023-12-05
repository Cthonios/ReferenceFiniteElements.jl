"""
"""
function element_stencil(::Tri6, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates = Ftype[
    0.0 1.0 0.0 0.5 0.5 0.0;
    0.0 0.0 1.0 0.0 0.5 0.5
  ]
  # @time face_nodes = Itype[
  #   1 2 3
  #   4 5 6
  #   2 3 1
  # ]
  face_nodes = Matrix{Itype}(undef, 3, 3)
  face_nodes[1, 1] = 1
  face_nodes[2, 1] = 4
  face_nodes[3, 1] = 2
  #
  face_nodes[1, 2] = 2
  face_nodes[2, 2] = 5
  face_nodes[3, 2] = 3
  #
  face_nodes[1, 1] = 3
  face_nodes[2, 1] = 6
  face_nodes[3, 1] = 1
  
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, face_nodes, interior_nodes
end

# using Tri6(1) as a template
for type in types_to_generate_interpolants(Tri6(Val(1)))
  """
  """
  @eval function shape_function_values(e::Tri6, ::Type{$(type[1])}, ξ::A) where A <: AbstractArray{<:Number, 1}
    λ = 1. - ξ[1] - ξ[2]
    N = $(type[1]){num_nodes(e), eltype(ξ)}(
      λ * (2. * λ - 1.),
      ξ[1] * (2. * ξ[1] - 1.),
      ξ[2] * (2. * ξ[2] - 1.),
      4. * ξ[1] * λ,
      4. * ξ[1] * ξ[2],
      4. * ξ[2] * λ
    )
  end

  """
  """
  @eval function shape_function_gradients(e::Tri6, ::Type{$(type[2])}, ξ::A) where A <: AbstractArray{<:Number}
    λ = 1. - ξ[1] - ξ[2]
    ∇N_ξ = $(type[2]){num_nodes(e), num_dimensions(e), eltype(ξ), num_nodes(e) * num_dimensions(e)}(
      -1. * (2. * λ - 1.) - 2. * λ,
       (2. * ξ[1] - 1.) + 2. * ξ[1],
       0.0,
       4. * λ - 4. * ξ[1],
       4. * ξ[2],
      -4. * ξ[2],
      #
      -1. * (2. * λ - 1.) - 2. * λ,
       0.0,
       (2. * ξ[2] - 1.) + 2. * ξ[2],
      -4. * ξ[1],
       4. * ξ[1],
       4. * λ - 4. * ξ[2]
    )
  end

  """
  """
  @eval function shape_function_hessians(e::Tri6, ::Type{$(type[3])}, ξ::A) where A <: AbstractArray{<:Number, 1}
    N, D = num_nodes(e), num_dimensions(e)
    λ = 1. - ξ[1] - ξ[2]
    ∇∇N_ξ = $(type[3]){Tuple{N, D, D}, eltype(ξ), 3, N * D * D}(
      4., 4., 0., -8., 0., 0.,
      4., 0., 0., -4., 4., -4.,
      4., 0., 0., -4., 4., -4.,
      4., 0., 4.,  0., 0., -8.
    )
  end
end

