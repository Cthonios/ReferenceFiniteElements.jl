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
end

"""
"""
function shape_function_gradients(::Tri3, ξ::SVector{2, <:Real})
  ∇N_ξ = (@SMatrix [
    -1. -1.;
     1.  0.;
     0.  1.
  ]) |> SMatrix{3, 2, eltype(ξ), 6}
end

"""
"""
function shape_function_hessians(::Tri3, ξ::SVector{2, <:Real})
  ∇∇N_ξ = zeros(SArray{Tuple{3, 2, 2}, eltype(ξ), 3, 12})
end
