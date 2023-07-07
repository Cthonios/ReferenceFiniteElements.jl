"""
"""
struct SimplexTri <: TriUnion
end

"""
"""
function ReferenceFEStencil(e::SimplexTri, degree::I, Itype::Type = Integer, Rtype::Type = Float64) where I <: Integer
  n_points = (degree + 1) * (degree + 2) / 2 |> Integer
  
  # to be consistent with optimism
  degree = degree + 1

  lobatto_points, _ = gausslobatto(degree)
  lobatto_points .= (1. .+ lobatto_points) ./ 2.

  # points = zeros(Float64, 2, n_points)
  points = Matrix{Float64}(undef, 2, n_points)
  point = 1
  for i in 1:degree
    for j in 1:degree + 1 - i
      k = degree - i - j + 2
      points[1, point] = (1. + 2. * lobatto_points[k] - lobatto_points[j] - lobatto_points[i]) / 3.
      points[2, point] = (1. + 2. * lobatto_points[j] - lobatto_points[i] - lobatto_points[k]) / 3.
      point = point + 1
    end
  end

  vertex_points = [1, degree, n_points]

  ii = 0:degree - 1 |> collect
  jj = cumsum(reverse(ii)) + ii
  kk = reverse(jj) - ii

  ii .= ii .+ 1
  jj .= jj .+ 1
  kk .= kk .+ 1

  face_points = hcat(ii, jj, kk)
  interior_nodes = findall(x -> x ∉ face_points, 1:n_points)
  return ReferenceFEStencil{Itype, Rtype, SimplexTri}(e, degree - 1, points, vertex_points, face_points, interior_nodes)
end

"""
"""
function ShapeFunctions(
  ::SimplexTri, 
  stencil::ReferenceFEStencil{Itype, Rtype}, 
  q_rule::Quadrature{Rtype}
) where {Itype <: Integer, Rtype <: Real}

  n_nodes = (stencil.degree + 1) * (stencil.degree + 2) ÷ 2

  A  = zeros(Rtype, n_nodes, size(stencil.coordinates, 2))
  Ax = zeros(Rtype, n_nodes, size(stencil.coordinates, 2))
  Ay = zeros(Rtype, n_nodes, size(stencil.coordinates, 2))

  nf  = zeros(Rtype, n_nodes, size(q_rule.ξs, 2))
  nfx = zeros(Rtype, n_nodes, size(q_rule.ξs, 2))
  nfy = zeros(Rtype, n_nodes, size(q_rule.ξs, 2))

  vander2d!(A, Ax, Ay, stencil.coordinates, stencil.degree)
  vander2d!(nf, nfx, nfy, q_rule.ξs, stencil.degree)

  Ns = A \ nf
  dshapes_x = A \ nfx
  dshapes_y = A \ nfy

  ∇N_ξs = Array{Rtype, 3}(undef, size(Ns, 1), 2, size(Ns, 2))
  ∇N_ξs[:, 1, :] .= dshapes_x
  ∇N_ξs[:, 2, :] .= dshapes_y

  # convert to static arrays
  Ns    = SVector{n_nodes, Rtype}.(eachcol(Ns))
  ∇N_ξs = SMatrix{n_nodes, 2, Rtype}.(eachslice(∇N_ξs, dims=3))

  return StructArray{ShapeFunctionPair{n_nodes, 2, Rtype}}((Ns, ∇N_ξs))
end

# internals
function pascal_triangle_monomials(degree::Itype) where Itype <: Integer
  pq = Matrix{Integer}(undef, 2, sum(1:degree + 1))
  range_begin = 1
  for i in 1:degree + 1
    monomial_indices = 1:i
    pq[2, range_begin:range_begin + i - 1] = monomial_indices
    pq[1, range_begin:range_begin + i - 1] = reverse(monomial_indices)
    range_begin = range_begin + i
  end
  pq
end

function map_from_tri_to_square_old(ξs::Matrix{<:Real})
  small = 1e-12
  index_singular = ξs[:, 2] .> 1. - small
  ξs[index_singular, 2] .= 1. - small

  ηs = zeros(typeof(ξs[1]), size(ξs))
  ηs[:, 1] .= 2. * (1. .+ ξs[:, 1]) ./ (1. .- ξs[:, 2]) .- 1.
  ηs[:, 2] .= ξs[:, 2]
  ηs[index_singular, 1] .= -1.
  ηs[index_singular, 2] .= 1.

  dηs = zeros(typeof(ξs[1]), size(ξs))
  dηs[:, 1] = 2. ./ (1 .- ξs[:, 2])
  dηs[:, 2] = 2. * (1. .+ ξs[:, 1]) ./ (1. .- ξs[:, 2]).^2
  return ηs, dηs
end

function map_from_tri_to_square(ξs_in::Matrix{T}) where T <: Real
  small = 1e-12::T
  ξs = copy(ξs_in)

  index_singular = 2
  ξs[2, index_singular] = 1. - small

  ηs = zeros(typeof(ξs[1]), size(ξs))
  # @time @views ηs[1, :] .= 2. * (1. .+ ξs[1, :]) ./ (1. .- ξs[2, :]) .- 1.
  @views @. ηs[1, :] = 2. * (1. + ξs[1, :]) / (1. - ξs[2, :]) - 1.
  @views ηs[2, :] .= ξs[2, :]

  ηs[1, index_singular] = -1.
  ηs[2, index_singular] = 1.

  dηs = zeros(typeof(ξs[1]), size(ξs))
  @views @. dηs[1, :] = 2. / (1 - ξs[2, :])
  @views @. dηs[2, :] = 2. * (1. + ξs[1, :]) / (1. - ξs[2, :])^2
  return ηs, dηs
end

function vander2d!(
  A::Matrix{Rtype}, Ax::Matrix{Rtype}, Ay::Matrix{Rtype},
  Xs::Matrix{Rtype}, degree::Itype
) where {Itype <: Integer, Rtype <: Real}
  n_nodes = (degree + 1) * (degree + 2) ÷ 2

  pq = pascal_triangle_monomials(degree)

  zs = @. 2. * Xs - 1.

  # ηs, dηs = map_from_tri_to_square(zs)
  ηs, dηs = map_from_tri_to_square_old(zs' |> collect)
  ηs = ηs'
  dηs = dηs'

  N1D = Polynomial([0.5, -0.5])
  for i in 1:n_nodes
    p = basis(Legendre, pq[1, i] - 1)
    temp = zeros(Rtype, pq[2, i])
    temp[end] = 1
    q = Jacobi{2 * (pq[1, i] - 1) + 1, 0}(temp)
    for _ in 1:pq[1, i] - 1
      q = q * N1D
    end
    dp = derivative(p)
    dq = derivative(q)

    @views weight = sqrt((2 * (pq[1, i] - 1) + 1) * 2 * ((pq[1, i] - 1) + (pq[2, i] - 1) + 1))

    @views A[i, :]  .= weight * p.(ηs[1, :]) .* q.(ηs[2, :])
    @views Ax[i, :] .= 2 * weight * dp.(ηs[1, :]) .* q.(ηs[2, :]) .* dηs[1, :]
    @views Ay[i, :] .= 2 * weight * (dp.(ηs[1, :]) .* q.(ηs[2, :]) .* dηs[2, :] + p.(ηs[1, :]) .* dq.(ηs[2, :]))
  end
end

export SimplexTri
