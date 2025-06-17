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

function CellInterpolants(e::AbstractElementType, Xs, backend::ArrayBackend, Ftype::Type{T}) where T <: Number
  ξs, ws = quadrature_points_and_weights(e, backend)
  Ns = map(x -> shape_function_value(e, Xs, x, backend), ξs)
  ∇N_ξs = map(x -> shape_function_gradient(e, Xs, x, backend), ξs)
  ∇∇N_ξs = map(x -> shape_function_hessian(e, Xs, x, backend), ξs)

  D = dimension(e)
  N = num_vertices(e)

  nt = NamedTuple{
    (:ξ, :w, :N, :∇N_ξ, :∇∇N_ξ),
    Tuple{
      Vector{SVector{D, Ftype}},
      Vector{Ftype},
      Vector{SVector{N, Ftype}},
      Vector{SMatrix{N, D, Ftype, N * D}},
      Vector{SArray{Tuple{N, D, D}, Ftype, 3, N * D * D}}
    }
  }((ξs, ws, Ns, ∇N_ξs, ∇∇N_ξs))
  vals = StructArray{Interpolants, 1, typeof(nt)}(nt)
  return CellInterpolants{typeof(vals)}(vals)
end

function SurfaceInterpolants(e::AbstractElementType, Xs, backend::ArrayBackend, Ftype::Type{T}) where T <: Number
  ξs, ws = surface_quadrature_points_and_weights(e, backend)
  ξs = mapreduce(x -> x, hcat, ξs)
  ws = mapreduce(x -> x, hcat, ws)
  Ns = shape_function_value.((e,), (Xs,), ξs, (backend,))
  ∇N_ξs = shape_function_gradient.((e,), (Xs,), ξs, (backend,))
  ∇∇N_ξs = shape_function_hessian.((e,), (Xs,), ξs, (backend,))

  D = dimension(e)
  N = num_vertices(e)

  nt = NamedTuple{
    (:ξ, :w, :N, :∇N_ξ, :∇∇N_ξ),
    Tuple{
      Matrix{SVector{D, Ftype}},
      Matrix{Ftype},
      Matrix{SVector{N, Ftype}},
      Matrix{SMatrix{N, D, Ftype, N * D}},
      Matrix{SArray{Tuple{N, D, D}, Ftype, 3, N * D * D}}
    }
  }((ξs, ws, Ns, ∇N_ξs, ∇∇N_ξs))
  vals = StructArray{Interpolants, 2, typeof(nt)}(nt)
  return SurfaceInterpolants{typeof(vals)}(vals)
end
