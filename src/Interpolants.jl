abstract type AbstractInterpolants end

# very general form
struct Interpolants{A, B, C, D, E} <: AbstractInterpolants
  ξ::A
  w::B
  N::C
  ∇N_ξ::D
  ∇∇N_ξ::E
end

# specific to staticarray backend
struct SInterpolants{RT <: Number, ND, NI, NIND, NINDND} <: AbstractInterpolants
  ξ::SVector{ND, RT}
  w::RT
  N::SVector{NI, RT}
  ∇N_ξ::SMatrix{NI, ND, RT, NIND}
  ∇∇N_ξ::SArray{Tuple{NI, ND, ND}, RT, 3, NINDND}
end

# special case for 1D elements.
function SInterpolants(
  ξ::RT, w::RT, 
  N::SVector{NI, RT}, 
  ∇N_ξ::SMatrix{NI, ND, RT, NIND},
  ∇∇N_ξ::SArray{Tuple{NI, ND, ND}, RT, 3, NINDND}
) where {RT, NI, ND, NIND, NINDND}
  return SInterpolants(SVector{ND, RT}(ξ), w, N, ∇N_ξ, ∇∇N_ξ)
end

function SInterpolants(
  ξ::AbstractVector{RT}, w::RT, 
  N::SVector{NI, RT}, 
  ∇N_ξ::SMatrix{NI, ND, RT, NIND},
  ∇∇N_ξ::SArray{Tuple{NI, ND, ND}, RT, 3, NINDND}
) where {RT, NI, ND, NIND, NINDND}
  return SInterpolants(SVector{ND, RT}(ξ...), w, N, ∇N_ξ, ∇∇N_ξ)
end

# function Base.show(io::IO, interps::AbstractInterpolantsContainer; tab="")
#   println(io, "$(tab)CellInterpolants Type                         = $(typeof(interps))")
#   # println(io, "$(tab)  Quadrature Point Storage Type               = $(typeof(interps.vals[1].ξ))")
#   # println(io, "$(tab)  Quadrature Weights Storage Type             = $(typeof(interps.vals[1].w))")
#   # println(io, "$(tab)  Shape Function Values Storage Type          = $(typeof(interps.vals[1].N))")
#   # println(io, "$(tab)  Shape Function Gradients Storage Type       = $(typeof(interps.vals[1].∇N_ξ))")
#   # println(io, "$(tab)  Shape Function Hessians Storage Type        = $(typeof(interps.vals[1].∇∇N_ξ))")
# end

# struct CellInterpolants{I} <: AbstractInterpolantsContainer{I}
#   vals::I
# end

# struct SurfaceInterpolants{I} <: AbstractInterpolantsContainer{I}
#   vals::I
# end

function _setup_cell_interpolants(e::AbstractElementType, Xs, backend::ArrayBackend)
  ξs, ws = quadrature_points_and_weights(e, backend)
  Ns = map(x -> shape_function_value(e, Xs, x, backend), ξs)
  ∇N_ξs = map(x -> shape_function_gradient(e, Xs, x, backend), ξs)
  ∇∇N_ξs = map(x -> shape_function_hessian(e, Xs, x, backend), ξs)

  if get_backend(backend) == SArray
    RT = eltype(ws)
    ND = dimension(e)
    NI = length(Ns[1])
    type = SInterpolants{RT, ND, NI, NI * ND, NI * ND * ND}
  else
    @assert false "Unsupported backend $backend"
  end

  vals = map((a, b, c, d, e) -> type(a, b, c, d, e), ξs, ws, Ns, ∇N_ξs, ∇∇N_ξs)

  return vals
end

function _setup_surface_interpolants(e::AbstractElementType, Xs, backend::ArrayBackend)
  ξs, ws = surface_quadrature_points_and_weights(e, backend)
  ξs = mapreduce(x -> x, hcat, ξs)
  ws = mapreduce(x -> x, hcat, ws)
  Ns = shape_function_value.((e,), (Xs,), ξs, (backend,))
  ∇N_ξs = shape_function_gradient.((e,), (Xs,), ξs, (backend,))
  ∇∇N_ξs = shape_function_hessian.((e,), (Xs,), ξs, (backend,))

  RT = eltype(ws)
  ND = dimension(e)
  NI = length(Ns[1])
  # type = SInterpolants{RT, ND, NI, NI * ND, NI * ND * ND}

  vals = map((a, b, c, d, e) -> SInterpolants(a, b, c, d, e), ξs, ws, Ns, ∇N_ξs, ∇∇N_ξs)
  return vals
end
