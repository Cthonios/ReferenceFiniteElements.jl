"""
"""
function element_stencil(::Tri3, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  nodal_coordinates = Ftype[
    0.0 1.0 0.0;
    0.0 0.0 1.0
  ]
  face_nodes = Itype[
    1 2 3
    2 3 1
  ]
  interior_nodes = Vector{Itype}(undef, 0)
  return nodal_coordinates, face_nodes, interior_nodes
end

# using Tri3(1) as a template
for type in types_to_generate_interpolants(Tri3(1))
  """
  """
  @eval function shape_function_values(e::Tri3, ::Type{$(type[1])}, ξ::A) where A <: AbstractArray{<:Number, 1}
    N = $(type[1]){num_nodes(e), eltype(ξ)}(
      1. - ξ[1] - ξ[2],
      ξ[1],
      ξ[2]
    )
  end

  """
  """
  @eval function shape_function_gradients(e::Tri3, ::Type{$(type[2])}, ξ::A) where A <: AbstractArray{<:Number, 1}
    ∇N_ξ = $(type[2]){num_nodes(e), num_dimensions(e), eltype(ξ), num_nodes(e) * num_dimensions(e)}(
      -1., 1., 0.,
      -1., 0., 1.
    )
  end

  """
  """
  @eval function shape_function_hessians(e::Tri3, ::Type{$(type[3])}, ξ::A) where A <: AbstractArray{<:Number, 1}
    N, D = num_nodes(e), num_dimensions(e)
    ∇∇N_ξ = zeros($(type[3]){Tuple{N, D, D}, eltype(ξ), 3, N * D * D})
  end
end


