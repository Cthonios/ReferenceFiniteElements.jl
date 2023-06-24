struct ShapeFunctionPair{N, D, Rtype <: Real}
  N::SVector{N, Rtype}
  ∇N_ξ::SMatrix{N, D, Rtype}
end

function ShapeFunctions end

function ShapeFunctions(
  e::E, 
  q_rule::Quadrature{Rtype}
) where {E <: ReferenceFE, Rtype <: Real}
  Ns = shape_function_values.((e,), eachcol(q_rule.ξs))
  ∇N_ξs = shape_function_gradients.((e,), eachcol(q_rule.ξs))
  n_nodes, n_dims = size(∇N_ξs[1])
  return StructArray{ShapeFunctionPair{n_nodes, n_dims, Rtype}}((Ns, ∇N_ξs))
end

function ShapeFunctions(e::E, degree::Int) where {E <: ReferenceFE}
  return ShapeFunctions(
    e, 
    # ElementStencil(e, degree),
    Quadrature(e, degree)
  )
end

export ShapeFunctions # this is really a method defined elsewhere
export ShapeFunctionPair
