function cell_quadrature_points_and_weights(::AbstractTri, q_rule::GaussLobattoLegendre)
    if cell_quadrature_degree(q_rule) == 1
      ξs = Matrix{Float64}(undef, 2, 1)
      ξs[:, 1] = [1. / 3., 1. / 3.]
      ws = [0.5]
    elseif cell_quadrature_degree(q_rule) == 2
      ξs = Matrix{Float64}(undef, 2, 3)
      ξs[:, 1] = [2. / 3., 1. / 6.]
      ξs[:, 2] = [1. / 6., 2. / 3.]
      ξs[:, 3] = [1. / 6., 1. / 6.]
      ws = [1. / 6., 1. / 6., 1. / 6.]
    elseif cell_quadrature_degree(q_rule) <= 4
      ξs = Matrix{Float64}(undef, 2, 6)
      ξs[:, 1] = [1.081030181680700E-01, 4.459484909159650E-01]
      ξs[:, 2] = [4.459484909159650E-01, 1.081030181680700E-01]
      ξs[:, 3] = [4.459484909159650E-01, 4.459484909159650E-01]
      ξs[:, 4] = [8.168475729804590E-01, 9.157621350977100E-02]
      ξs[:, 5] = [9.157621350977100E-02, 8.168475729804590E-01]
      ξs[:, 6] = [9.157621350977100E-02, 9.157621350977100E-02]
  
      ws = [
        1.116907948390055E-01,
        1.116907948390055E-01,
        1.116907948390055E-01,
        5.497587182766100E-02,
        5.497587182766100E-02,
        5.497587182766100E-02
      ]
    elseif cell_quadrature_degree(q_rule) <= 5
      ξs = Matrix{Float64}(undef, 2, 7)
      ξs[:, 1] = [3.33333333333333E-01, 3.33333333333333E-01]
      ξs[:, 2] = [5.97158717897700E-02, 4.70142064105115E-01]
      ξs[:, 3] = [4.70142064105115E-01, 5.97158717897700E-02]
      ξs[:, 4] = [4.70142064105115E-01, 4.70142064105115E-01]
      ξs[:, 5] = [7.97426985353087E-01, 1.01286507323456E-01]
      ξs[:, 6] = [1.01286507323456E-01, 7.97426985353087E-01]
      ξs[:, 7] = [1.01286507323456E-01, 1.01286507323456E-01]
  
      ws = [
        1.12500000000000E-01,
        6.61970763942530E-02,
        6.61970763942530E-02,
        6.61970763942530E-02,
        6.29695902724135E-02,
        6.29695902724135E-02,
        6.29695902724135E-02
      ]
    elseif cell_quadrature_degree(q_rule) <= 6
      ξs = Matrix{Float64}(undef, 2, 12)
      ξs[:, 1]  = [5.01426509658179E-01, 2.49286745170910E-01]
      ξs[:, 2]  = [2.49286745170910E-01, 5.01426509658179E-01]
      ξs[:, 3]  = [2.49286745170910E-01, 2.49286745170910E-01]
      ξs[:, 4]  = [8.73821971016996E-01, 6.30890144915020E-02]
      ξs[:, 5]  = [6.30890144915020E-02, 8.73821971016996E-01]
      ξs[:, 6]  = [6.30890144915020E-02, 6.30890144915020E-02]
      ξs[:, 7]  = [5.31450498448170E-02, 3.10352451033784E-01]
      ξs[:, 8]  = [6.36502499121399E-01, 5.31450498448170E-02]
      ξs[:, 9]  = [3.10352451033784E-01, 6.36502499121399E-01]
      ξs[:, 10] = [5.31450498448170E-02, 6.36502499121399E-01]
      ξs[:, 11] = [6.36502499121399E-01, 3.10352451033784E-01]
      ξs[:, 12] = [3.10352451033784E-01, 5.31450498448170E-02]
  
      ws = [
        5.83931378631895E-02,
        5.83931378631895E-02,
        5.83931378631895E-02,
        2.54224531851035E-02,
        2.54224531851035E-02,
        2.54224531851035E-02,
        4.14255378091870E-02,
        4.14255378091870E-02,
        4.14255378091870E-02,
        4.14255378091870E-02,
        4.14255378091870E-02,
        4.14255378091870E-02
      ]
    else
        @assert false "Quadrature degree 1 through 6 currently supported."
    end

    return ξs, ws
end

function surface_quadrature_points_and_weights(e::AbstractTri, q_rule::GaussLobattoLegendre)
    ξs, ws = cell_quadrature_points_and_weights(boundary_element(e, 0), q_rule)

    ξ_return = zeros(2, length(ws), 3)
    w_return = zeros(length(ws), 3)

    ξ_return[1, :, 1] .= ξs[1, :]
    ξ_return[2, :, 1] .= -1.
    ξ_return[1, :, 2] .= ξs[1, :]
    ξ_return[2, :, 2] .= 1. .- ξs[1, :]
    ξ_return[1, :, 3] .= -1.
    ξ_return[2, :, 3] .= ξs[1, :]

    for n in 1:3
        w_return[:, n] .= ws
    end
    return ξ_return, w_return
end

struct Tri{PT, PD} <: AbstractTri{PT, PD}
end

########################################################################
# Lagrange implementation
########################################################################
function boundary_dofs(e::Tri{Lagrange, PD}) where PD
    linear_edges = edge_vertices(e)
    if PD < 2
        return linear_edges
    else
        edges = zeros(Int, PD + 1, 3)
        edges[1:2, 1:3] .= linear_edges
        offset = 4
        for n in 1:3
            edges[3:end, n] = offset:offset + PD - 2
            offset += PD - 1
        end
        return edges
    end
end
function dof_coordinates(e::Tri{Lagrange, PD}) where PD
    if PD == 0
        return zeros(2, 1)
    end

    coords = [
        0. 1. 0.;
        0. 0. 1.
    ]

    if PD > 1
        # do edge midpoints
        edge_coords = dof_coordinates(boundary_element(e, 0))

        # face 1
        for n in 1:PD - 1
            coords = hcat(coords, [edge_coords[1, n + 2], -1.])
        end
        # face 2
        for n in 1:PD - 1
            coords = hcat(coords, [edge_coords[1, n + 2], 1. - edge_coords[1, n + 2]])
        end
        # face 3
        for n in 1:PD - 1
            coords = hcat(coords, [-1., edge_coords[1, n + 2]])
        end

        # now for interiors
        # TODO not correct yet
        for n in 1:PD - 1
            for m in 1:PD - 1 - n
                # @assert false "TODO"
                coords = hcat(coords, [edge_coords[1, m + 2], edge_coords[1, n + 2]])
            end
        end
    end
    return coords
end
function interior_dofs(::Tri{Lagrange, PD}) where PD
    if PD < 3
        return Int[]
    else
        @assert false "TODO"
    end
end
num_cell_dofs(::Tri{Lagrange, PD}) where PD = (PD + 1) * (PD + 2) ÷ 2
num_interior_dofs(::Tri{Lagrange, PD}) where PD = PD < 3 ? 0 : (PD - 1) * (PD - 2) ÷ 2

function shape_function_value(::Tri{Lagrange, 0}, _, _)
    return ones(1)
end

function shape_function_gradient(::Tri{Lagrange, 0}, _, _)
    return zeros(1, 2)
end

function shape_function_hessian(::Tri{Lagrange, 0}, _, _)
    return zeros(1, 2, 2)
end

function shape_function_value(::Tri{Lagrange, PD}, _, ξ) where PD
    λ1 = 1 - ξ[1] - ξ[2]
    λ2 = ξ[1]
    λ3 = ξ[2]

    N = Vector{eltype(ξ)}(undef, (PD+1)*(PD+2)÷2)
    offset = 0

    # -------------------------
    # vertices
    # -------------------------
    offset += 1; N[offset] = λ1^PD
    offset += 1; N[offset] = λ2^PD
    offset += 1; N[offset] = λ3^PD

    # -------------------------
    # edges
    # -------------------------
    # edge 1–2 (λ3 = 0)
    for i in 1:PD - 1
        offset += 1
        N[offset] =
            binomial(PD, i) * λ1^(PD - i) * λ2^i
    end

    # edge 2–3 (λ1 = 0)
    for i in 1:PD - 1
        offset += 1
        N[offset] =
            binomial(PD, i) * λ2^(PD - i) * λ3^i
    end

    # edge 3–1 (λ2 = 0)
    for i in 1:PD - 1
        offset += 1
        N[offset] =
            binomial(PD, i) * λ3^(PD - i) * λ1^i
    end

    # -------------------------
    # interior
    # -------------------------
    for i in 1:PD - 2, j in 1:PD - 1 - i
        k = PD - i - j
        offset += 1
        N[offset] =
            binomial(PD, i) *
            binomial(PD - i, j) *
            λ1^i * λ2^j * λ3^k
    end

    return N
end

function shape_function_gradient(e::Tri{Lagrange, PD}, _, ξ) where PD
    λ1 = 1 - ξ[1] - ξ[2]
    λ2 = ξ[1]
    λ3 = ξ[2]

    gλ1 = SVector(-1.0, -1.0)
    gλ2 = SVector(1.0,  0.0)
    gλ3 = SVector(0.0,  1.0)

    ndofs = num_cell_dofs(e)
    dN = Matrix{eltype(ξ)}(undef, ndofs, 2)
    offset = 0

    # -------------------------
    # vertices
    # -------------------------
    offset += 1
    dN[offset, :] = PD * λ1^(PD - 1) * gλ1

    offset += 1
    dN[offset, :] = PD * λ2^(PD - 1) * gλ2

    offset += 1
    dN[offset, :] = PD * λ3^(PD - 1) * gλ3

    # -------------------------
    # edges
    # -------------------------
    # edge 1–2 (λ3 = 0)
    for i in 1:PD - 1
        offset += 1
        C = binomial(PD, i)
        dN[offset, :] =
            C * (
                (PD - i) * λ1^(PD - i - 1) * λ2^i       * gλ1 +
                i        * λ1^(PD - i)     * λ2^(i - 1) * gλ2
            )
    end

    # edge 2–3 (λ1 = 0)
    for i in 1:PD - 1
        offset += 1
        C = binomial(PD, i)
        dN[offset, :] =
            C * (
                (PD - i) * λ2^(PD - i - 1) * λ3^i       * gλ2 +
                i        * λ2^(PD - i)     * λ3^(i - 1) * gλ3
            )
    end

    # edge 3–1 (λ2 = 0)
    for i in 1:PD - 1
        offset += 1
        C = binomial(PD, i)
        dN[offset, :] =
            C * (
                (PD - i) * λ3^(PD - i - 1) * λ1^i       * gλ3 +
                i        * λ3^(PD - i)     * λ1^(i - 1) * gλ1
            )
    end

    # -------------------------
    # interior
    # -------------------------
    for i in 1:PD - 2, j in 1:PD - 1 - i
        k = PD - i - j
        offset += 1

        C = binomial(PD, i) * binomial(PD - i, j)

        dN[offset, :] =
            C * (
                i * λ1^(i - 1) * λ2^j       * λ3^k       * gλ1 +
                j * λ1^i       * λ2^(j - 1) * λ3^k       * gλ2 +
                k * λ1^i       * λ2^j       * λ3^(k - 1) * gλ3
            )
    end

    return dN
end

function shape_function_hessian(e::Tri{Lagrange, PD}, _, ξ) where PD
    T = eltype(ξ)

    λ1 = one(T) - ξ[1] - ξ[2]
    λ2 = ξ[1]
    λ3 = ξ[2]

    gλ1 = SVector(-one(T), -one(T))
    gλ2 = SVector(one(T), zero(T))
    gλ3 = SVector(zero(T), one(T))

    ndofs = num_cell_dofs(e)
    H = Array{T}(undef, ndofs, 2, 2)

    @inline pow0(λ, p) = p == 0 ? one(T) : p > 0 ? λ^p : zero(T)

    @inline sym(a, b) = SMatrix{2,2}(
        a[1] * b[1], a[1] * b[2],
        a[2] * b[1], a[2] * b[2]
    )

    offset = 0
    for i in 0:PD, j in 0:(PD - i)
        k = PD - i - j
        offset += 1
        H[offset, :, :] .= zero(T)

        C = binomial(PD, i) * binomial(PD - i, j)

        if i ≥ 2
            H[offset, :, :] .+=
                C * i * (i - 1) *
                pow0(λ1, i - 2) * pow0(λ2, j) * pow0(λ3, k) * sym(gλ1, gλ1)
        end
        if j ≥ 2
            H[offset, :, :] .+=
                C * j * (j - 1) *
                pow0(λ1, i) * pow0(λ2, j - 2) * pow0(λ3, k) * sym(gλ2, gλ2)
        end
        if k ≥ 2
            H[offset, :, :] .+=
                C * k * (k - 1) *
                pow0(λ1, i) * pow0(λ2, j) * pow0(λ3, k - 2) * sym(gλ3, gλ3)
        end

        if i ≥ 1 && j ≥ 1
            H[offset, :, :] .+=
                C * i * j *
                pow0(λ1, i - 1) * pow0(λ2, j - 1) * pow0(λ3, k) * (sym(gλ1, gλ2) + sym(gλ2, gλ1))
        end
        if i ≥ 1 && k ≥ 1
            H[offset, :, :] .+=
                C * i * k *
                pow0(λ1, i-1) * pow0(λ2, j) * pow0(λ3, k - 1) * (sym(gλ1, gλ3) + sym(gλ3, gλ1))
        end
        if j ≥ 1 && k ≥ 1
            H[offset, :, :] .+=
                C * j * k *
                pow0(λ1, i) * pow0(λ2, j - 1) * pow0(λ3, k - 1) * (sym(gλ2, gλ3) + sym(gλ3, gλ2))
        end
    end

    return H
end

########################################################################
# Raviart-Thomas implementation
########################################################################
function boundary_dofs(::Tri{RaviartThomas, 0})
    return reshape(collect(1:3), 1, 3)
end
function dof_coordinates(::Tri{RaviartThomas, 0})
    # edge midpoints of reference triangle
    return [
        0.5  0.5  0.0;
        0.0  0.5  0.5
    ]
end
interior_dofs(::Tri{RaviartThomas, 0}) = Int[]
num_cell_dofs(::Tri{RaviartThomas, 0}) = 3
num_interior_dofs(::Tri{RaviartThomas, 0}) = 0

function geometry_shape_function_value(::Tri{RaviartThomas, 0}, X, ξ)
    return shape_function_value(Tri{Lagrange, 1}(), X, ξ)
end

function geometry_shape_function_gradient(::Tri{RaviartThomas, 0}, X, ξ)
    return shape_function_gradient(Tri{Lagrange, 1}(), X, ξ)
end

function shape_function_value(::Tri{RaviartThomas, 0}, _, ξ)
    N = Matrix{Float64}(undef, 3, 2)

    N[1, 1] = 1
    N[1, 2] = -ξ[2]
    #
    N[2, 1] = -ξ[1]
    N[2, 2] = 1 - ξ[2]
    #
    N[3, 1] = ξ[1]
    N[3, 2] = ξ[2]

    return N
end

function shape_function_divergence(::Tri{RaviartThomas, 0}, _, ξ)
    return [-2., -2., 2.]
end
