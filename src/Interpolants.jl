struct Interpolants{A, B, C, D, E}
  ξ::A
  w::B
  N::C
  ∇N_ξ::D
  ∇∇N_ξ::E
end

function Base.show(io::IO, interps::AbstractInterpolantsContainer; tab="")
  println(io, "$(tab)CellInterpolants:")
  println(io, "$(tab)  Quadrature Point Storage Type               = $(typeof(interps.vals.ξ))")
  println(io, "$(tab)  Quadrature Weights Storage Type             = $(typeof(interps.vals.w))")
  println(io, "$(tab)  Shape Function Values Storage Type          = $(typeof(interps.vals.N))")
  println(io, "$(tab)  Shape Function Gradients Storage Type       = $(typeof(interps.vals.∇N_ξ))")
  println(io, "$(tab)  Shape Function Hessians Storage Type        = $(typeof(interps.vals.∇∇N_ξ))")
end

struct CellInterpolants{I} <: AbstractInterpolantsContainer{I}
  vals::I
end

struct SurfaceInterpolants{I} <: AbstractInterpolantsContainer{I}
  vals::I
end

function CellInterpolants(e::AbstractElementType, Xs, backend::ArrayBackend)
  ξs, ws = quadrature_points_and_weights(e, backend)
  Ns = map(x -> shape_function_value(e, Xs, x, backend), ξs)
  ∇N_ξs = map(x -> shape_function_gradient(e, Xs, x, backend), ξs)
  ∇∇N_ξs = map(x -> shape_function_hessian(e, Xs, x, backend), ξs)
  vals = StructArray{Interpolants}((ξs, ws, Ns, ∇N_ξs, ∇∇N_ξs))
  return CellInterpolants(vals)
end

function SurfaceInterpolants(e::AbstractElementType, Xs, backend::ArrayBackend)
  ξs, ws = surface_quadrature_points_and_weights(e, backend)
  ξs = mapreduce(x -> x, hcat, ξs)
  ws = mapreduce(x -> x, hcat, ws)
  Ns = shape_function_value.((e,), (Xs,), ξs, (backend,))
  ∇N_ξs = shape_function_gradient.((e,), (Xs,), ξs, (backend,))
  ∇∇N_ξs = shape_function_hessian.((e,), (Xs,), ξs, (backend,))
  vals = StructArray{Interpolants}((ξs, ws, Ns, ∇N_ξs, ∇∇N_ξs))
  return CellInterpolants(vals)
end
