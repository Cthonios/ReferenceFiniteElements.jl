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

function shape_function_values(::Tri6, ξ::SVector{2, <:Real})
  λ = 1. - ξ[1] - ξ[2]
  N = SVector{6, eltype(ξ)}(
    λ * (2. * λ - 1.),
    ξ[1] * (2. * ξ[1] - 1.),
    ξ[2] * (2. * ξ[2] - 1.),
    4. * ξ[1] * λ,
    4. * ξ[1] * ξ[2],
    4. * ξ[2] * λ
  )
end

function shape_function_gradients(::Tri6, ξ::SVector{2, <:Real})
  λ = 1. - ξ[1] - ξ[2]
  ∇N_ξ = (@SMatrix [
    -1. * (2. * λ - 1.) - 2. * λ  -1. * (2. * λ - 1.) - 2. * λ;
     (2. * ξ[1] - 1.) + 2. * ξ[1]  0.0;
     0.0                           (2. * ξ[2] - 1.) + 2. * ξ[2];
     4. * λ - 4. * ξ[1]           -4. * ξ[1];
     4. * ξ[2]                     4. * ξ[1];
    -4. * ξ[2]                     4. * λ - 4. * ξ[2]
  ]) |> SMatrix{6, 2, eltype(ξ), 12}
end 

function shape_function_hessians(::Tri6, ξ::SVector{2, <:Real})
  λ = 1. - ξ[1] - ξ[2]
  ∇∇N_ξ = (@SArray [
    4.; 4.; 0.; -8.; 0.; 0.;;
    4.; 0.; 0.; -4.; 4.; -4.;;;
    4.; 0.; 0.; -4.; 4.; -4.;;
    4.; 0.; 4.;  0.; 0.; -8.;;;
  ]) |> SArray{Tuple{6, 2, 2}, eltype(ξ), 3, 24}
end