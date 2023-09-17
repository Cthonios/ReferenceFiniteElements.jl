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

"""
"""
function shape_function_values(::Quad4, ξ::SVector{2, <:Real})
  N = SVector{4, eltype(ξ)}(
    0.25 * (1.0 - ξ[1]) * (1.0 - ξ[2]),
    0.25 * (1.0 + ξ[1]) * (1.0 - ξ[2]),
    0.25 * (1.0 + ξ[1]) * (1.0 + ξ[2]),
    0.25 * (1.0 - ξ[1]) * (1.0 + ξ[2])
  )
end

"""
"""
function shape_function_gradients(::Quad4, ξ::SVector{2, <:Real})
  ∇N_ξ = (@SMatrix [
    -0.25 * (1.0 - ξ[2]) -0.25 * (1.0 - ξ[1]);
     0.25 * (1.0 - ξ[2]) -0.25 * (1.0 + ξ[1]);
     0.25 * (1.0 + ξ[2])  0.25 * (1.0 + ξ[1]);
    -0.25 * (1.0 + ξ[2])  0.25 * (1.0 - ξ[1])
  ]) |> SMatrix{4, 2, eltype(ξ), 8}
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
