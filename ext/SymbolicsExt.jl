module SymbolicsExt

using ReferenceFiniteElements
using Symbolics

function ReferenceFiniteElements.shape_function_value(e::Edge{Lagrange, PD}, X, ξ::Num) where PD
  Ns = Vector{Num}(undef, num_cell_dofs(e))
  for n in axes(Ns, 1)
    Ns[n] = 1
    for i in axes(Ns, 1)
      if i == n
        continue
      end
      Ns[n] = Ns[n] * (ξ[1] - X[1, i]) / (X[1, n] - X[1, i])
    end
  end

  return Ns
end

function ReferenceFiniteElements.shape_function_gradient(e::Edge, X, ξ::Num)
  Ns = ReferenceFiniteElements.shape_function_value(e, X, ξ)
  ∇N_ξs = Symbolics.derivative(Ns, ξ)
  return reshape(∇N_ξs, 1, length(∇N_ξs))
end

function ReferenceFiniteElements.shape_function_hessian(e::Edge, X, ξ::Num)
  Ns = ReferenceFiniteElements.shape_function_value(e, X, ξ)
  ∇∇N_ξs = Symbolics.derivative(Symbolics.derivative(Ns, ξ), ξ)
  return reshape(∇∇N_ξs, 1, 1, length(∇∇N_ξs))
end

end # module