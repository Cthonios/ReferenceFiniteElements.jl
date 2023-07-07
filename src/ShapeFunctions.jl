"""
"""
struct ShapeFunctions{N, D, T}
  Ns::Vector{SVector{N, T}}
  ∇N_ξs::Vector{SMatrix{N, D, T}}
end
"""
"""
shape_function_values(s::ShapeFunctions) = getfield(s, :Ns)
"""
"""
shape_function_gradient(s::ShapeFunctions, i::Integer) = getfield(s, :∇N_ξs)[i]
"""
"""
shape_function_gradients(s::ShapeFunctions) = getfield(s, :∇N_ξs)

"""
"""
function ShapeFunctions(
  e::ReferenceFEType{N, D}, 
  q_rule::Quadrature{Rtype}
) where {N, D, Rtype <: Real}

  Ns = Vector{SVector{N, Rtype}}(undef, num_q_points(q_rule))
  ∇N_ξs = Vector{SMatrix{N, D, Rtype}}(undef, num_q_points(q_rule))
  for (n, ξ) in enumerate(eachcol(quadrature_points(q_rule)))
    Ns[n] = shape_function_values_int(e, ξ)
    ∇N_ξs[n] = shape_function_gradients_int(e, ξ)
  end
  return ShapeFunctions{N, D, Rtype}(Ns, ∇N_ξs)
end

function ShapeFunctions(e::E, degree::Int) where {E <: ReferenceFEType}
  return ShapeFunctions(
    e, 
    Quadrature(e, degree)
  )
end
