"""
"""
function element_stencil(e::SimplexTri, ::Type{Itype}, ::Type{Ftype}) where {Itype <: Integer, Ftype <: AbstractFloat}
  n_points = (e.degree + 1) * (e.degree + 2) ÷ 2

  # to be consistent with optimism
  degree = e.degree + 1

  lobatto_points, _ = gausslobatto(degree)
  lobatto_points .= (1. .+ lobatto_points) ./ 2.

  nodal_coordinates = Matrix{Ftype}(undef, 2, n_points)
  point = 1
  for i in 1:degree
    for j in 1:degree + 1 - i
      k = degree - i - j + 2
      nodal_coordinates[1, point] = (1. + 2. * lobatto_points[k] - lobatto_points[j] - lobatto_points[i]) / 3.
      nodal_coordinates[2, point] = (1. + 2. * lobatto_points[j] - lobatto_points[i] - lobatto_points[k]) / 3.
      point = point + 1
    end
  end

  # TODO what about thsi guy?
  # vertex_points = [1, degree, n_points]

  # ii = 0:degree - 1 |> collect
  # jj = cumsum(reverse(ii)) + ii
  # kk = reverse(jj) - ii

  # ii .= ii .+ 1
  # jj .= jj .+ 1
  # kk .= kk .+ 1

  ii = 1:degree
  jj = cumsum(reverse(0:degree - 1)) + ii .+ 1
  kk = reverse(jj) - ii .+ 1

  edge_nodes = hcat(ii, jj, kk)
  edge_nodes = convert.(Itype, edge_nodes)
  face_nodes = Itype[;;]
  interior_nodes = findall(x -> x ∉ edge_nodes, 1:n_points)
  interior_nodes = convert.(Itype, interior_nodes)
  return nodal_coordinates, edge_nodes, face_nodes, interior_nodes
end

"""
"""
function setup_interpolants_simplex_tri!(Ns, ∇N_ξs, ∇∇N_ξs, e, ξs, coordinates)
  n_nodes = (e.degree + 1) * (e.degree + 2) ÷ 2
  Ftype   = eltype(coordinates)

  A  = zeros(Ftype, n_nodes, size(coordinates, 2))
  Ax = zeros(Ftype, n_nodes, size(coordinates, 2))
  Ay = zeros(Ftype, n_nodes, size(coordinates, 2))

  nf  = zeros(Ftype, n_nodes, size(ξs, 1))
  nfx = zeros(Ftype, n_nodes, size(ξs, 1))
  nfy = zeros(Ftype, n_nodes, size(ξs, 1))

  temp_ξs = mapreduce(permutedims, vcat, ξs)'

  vander2d!(A, Ax, Ay, coordinates, e.degree)
  vander2d!(nf, nfx, nfy, temp_ξs, e.degree)

  shapes    = A \ nf
  dshapes_x = A \ nfx
  dshapes_y = A \ nfy

  for n in 1:num_q_points(e)
    if eltype(Ns) <: Union{SVector, MVector}
      Ns[n] = eltype(Ns)(@views shapes[:, n])
      ∇N_ξs[n] = eltype(∇N_ξs)(@views hcat(dshapes_x[:, n], dshapes_y[:, n]))
    else
      Ns[n] = shapes[:, n]
      ∇N_ξs[n] = hcat(dshapes_x[:, n], dshapes_y[:, n])'
    end
  end
end

function vander2d!(
  A::Matrix{Ftype}, Ax::Matrix{Ftype}, Ay::Matrix{Ftype},
  Xs, degree::Itype
) where {Itype <: Integer, Ftype <: Real}

  n_nodes = (degree + 1) * (degree + 2) ÷ 2
  pq = pascal_triangle_monomials(degree)

  zs = 2. * Xs .- 1.

  # ηs, dηs = map_from_tri_to_square(zs)
  ηs, dηs = map_from_tri_to_square(zs' |> collect)
  ηs = ηs'
  dηs = dηs'

  N1D = Polynomial([0.5, -0.5])
  for i in 1:n_nodes
    p = basis(Legendre, pq[1, i] - 1)
    temp = zeros(Ftype, pq[2, i])
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

function map_from_tri_to_square(ξs::Matrix{<:Real})
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