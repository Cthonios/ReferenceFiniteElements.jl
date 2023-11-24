"""
"""
function element_stencil(::Quad4, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
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

# using Quad4(1) as a template
for type in types_to_generate_interpolants(Quad4(1))
  """
  """
  @eval function shape_function_values(e::Quad4, ::Type{$(type[1])}, ξ::A) where A <: AbstractArray{<:Number, 1}
    N = $(type[1]){num_nodes(e), eltype(ξ)}(
      0.25 * (1.0 - ξ[1]) * (1.0 - ξ[2]),
      0.25 * (1.0 + ξ[1]) * (1.0 - ξ[2]),
      0.25 * (1.0 + ξ[1]) * (1.0 + ξ[2]),
      0.25 * (1.0 - ξ[1]) * (1.0 + ξ[2])
    )
  end

  """
  """
  @eval function shape_function_gradients(e::Quad4, ::Type{$(type[2])}, ξ::A) where A <: AbstractArray{<:Number, 1}
    ∇N_ξ = $(type[2]){num_nodes(e), num_dimensions(e), eltype(ξ), num_nodes(e) * num_dimensions(e)}(
      -0.25 * (1.0 - ξ[2]),
       0.25 * (1.0 - ξ[2]),
       0.25 * (1.0 + ξ[2]),
      -0.25 * (1.0 + ξ[2]),
      #
      -0.25 * (1.0 - ξ[1]),
      -0.25 * (1.0 + ξ[1]),
       0.25 * (1.0 + ξ[1]),
       0.25 * (1.0 - ξ[1])
    )
  end
end

"""
"""
function shape_function_hessians(::Quad4, ξ::SVector{2, <:Real})
  ∇∇N_ξ = (@SArray [
    0.0; 0.0; 0.0; 0.0;; 
    0.25; -0.25; 0.25; -0.25;;;
    0.25; -0.25; 0.25; -0.25;; 
    0.0; 0.0; 0.0; 0.0
  ]) |> SArray{Tuple{4, 2, 2}, eltype(ξ), 3, 16}
end
