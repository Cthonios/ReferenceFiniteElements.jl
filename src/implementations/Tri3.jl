"""
"""
function element_stencil(::Tri3, degree::I, ::Type{Itype}, ::Type{Ftype}) where {I <: Integer, Itype <: Integer, Ftype <: AbstractFloat}
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

function shape_function_values(::Tri3, ξ::T) where T <: AbstractArray
  N = @SVector [
    1. - ξ[1] - ξ[2],
    ξ[1],
    ξ[2]
  ]
end

function shape_function_gradients(::Tri3, ξ::T) where T <: AbstractArray
  ∇N_ξ = @SMatrix [
    -1. -1.;
     1.  0.;
     0.  1.
  ]
end
