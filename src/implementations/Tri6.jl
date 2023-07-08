"""
"""
function element_stencil(::Tri6, degree::I, ::Type{Itype}, ::Type{Ftype}) where {I <: Integer, Itype <: Integer, Ftype <: AbstractFloat}
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

function shape_function_values(::Tri6, ξ::T) where T <: AbstractArray
  λ = 1. - ξ[1] - ξ[2]
  N = @SVector [
    λ * (2. * λ - 1.),
    ξ[1] * (2. * ξ[1] - 1.),
    ξ[2] * (2. * ξ[2] - 1.),
    4. * ξ[1] * λ,
    4. * ξ[1] * ξ[2],
    4. * ξ[2] * λ
  ]
end

function shape_function_gradients(::Tri6, ξ::T) where T <: AbstractArray
  λ = 1. - ξ[1] - ξ[2]
  ∇N_ξ = @SMatrix [
    -1. * (2. * λ - 1.) - 2. * λ       -1. * (2. * λ - 1.) - 2. * λ      ;
     (2. * ξ[1] - 1.) + 2. * ξ[1]       0.0                              ;
     0.0                                (2. * ξ[2] - 1.) + 2. * ξ[2]     ;
     4. * λ - 4. * ξ[1]                -4. * ξ[1]                        ;
     4. * ξ[2]                          4. * ξ[1]                        ;
    -4. * ξ[2]                          4. * λ - 4. * ξ[2]               ;  
  ]
end 
