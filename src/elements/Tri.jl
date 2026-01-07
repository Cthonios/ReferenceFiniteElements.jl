struct Tri{PT, PD} <: AbstractTri{PT, PD}
end

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
        edge_coords = dof_coordinates(boundary_element(e))

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
    ξs, ws = cell_quadrature_points_and_weights(boundary_element(e), q_rule)

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

# TODO eventually use edge methods instead of the below
function _simplex_lagrange_1d(i::Int, λ, PD::Int)
    i == 0 && return one(λ)
    val = one(λ)
    xi = i / PD
    for m in 0:(i - 1)
        xm = m / PD
        val *= (λ - xm) / (xi - xm)
    end
    return val
end

function _simplex_lagrange_1d_derivative(i::Int, λ, PD::Int)
    i == 0 && return zero(λ)

    xi = i / PD
    dval = zero(λ)

    for n in 0:(i - 1)
        xn = n / PD

        term = one(λ)
        for m in 0:(i - 1)
            m == n && continue
            xm = m / PD
            term *= (λ - xm) / (xi - xm)
        end

        dval += term / (xi - xn)
    end

    return dval
end

function _simplex_lagrange_1d_second_derivative(i::Int, λ, PD::Int)
    i ≤ 1 && return zero(λ)

    xi = i / PD
    d2 = zero(λ)

    for a in 0:(i - 1), b in 0:(i - 1)
        a == b && continue

        xa = a / PD
        xb = b / PD

        term = one(λ)
        for m in 0:(i - 1)
            (m == a || m == b) && continue
            xm = m / PD
            term *= (λ - xm) / (xi - xm)
        end

        d2 += term / ((xi - xa) * (xi - xb))
    end

    return d2
end

function shape_function_value(e::Tri{Lagrange, PD}, _, ξ) where PD
    λ1 = 1.0 - ξ[1] - ξ[2]
    λ2 = ξ[1]
    λ3 = ξ[2]

    N = Vector{eltype(ξ)}(undef, num_cell_dofs(e))
    offset = 0

    # -------------------------
    # vertices
    # -------------------------
    N[1] = _simplex_lagrange_1d(PD, λ1, PD)
    N[2] = _simplex_lagrange_1d(PD, λ2, PD)
    N[3] = _simplex_lagrange_1d(PD, λ3, PD)
    offset = 3

    # -------------------------
    # edges
    # -------------------------
    # edge 1–2 (λ3 = 0)
    for i in 1:PD - 1
        offset += 1
        N[offset] =
            _simplex_lagrange_1d(PD - i, λ1, PD) *
            _simplex_lagrange_1d(i,      λ2, PD)
    end

    # edge 2–3 (λ1 = 0)
    for i in 1:PD - 1
        offset += 1
        N[offset] =
            _simplex_lagrange_1d(PD - i, λ2, PD) *
            _simplex_lagrange_1d(i,      λ3, PD)
    end

    # edge 3–1 (λ2 = 0)
    for i in 1:PD - 1
        offset += 1
        N[offset] =
            _simplex_lagrange_1d(PD - i, λ3, PD) *
            _simplex_lagrange_1d(i,      λ1, PD)
    end

    # -------------------------
    # interior
    # -------------------------
    for i in 1:PD - 2, j in 1:PD - 1 - i
        k = PD - i - j
        offset += 1
        N[offset] =
            _simplex_lagrange_1d(i, λ1, PD) *
            _simplex_lagrange_1d(j, λ2, PD) *
            _simplex_lagrange_1d(k, λ3, PD)
    end

    return N
end

function shape_function_gradient(e::Tri{Lagrange, PD}, _, ξ) where PD
    λ1 = 1.0 - ξ[1] - ξ[2]
    λ2 = ξ[1]
    λ3 = ξ[2]

    dN = Matrix{eltype(ξ)}(undef, 2, num_cell_dofs(e))
    offset = 0

    # barycentric gradients
    gλ1 = [-1., -1.]
    gλ2 = [1., 0.]
    gλ3 = [0., 1.]

    # -------------------------
    # vertices
    # -------------------------
    dN[:, 1] = _simplex_lagrange_1d_derivative(PD, λ1, PD) * gλ1
    dN[:, 2] = _simplex_lagrange_1d_derivative(PD, λ2, PD) * gλ2
    dN[:, 3] = _simplex_lagrange_1d_derivative(PD, λ3, PD) * gλ3

    offset = 3

    # -------------------------
    # edges
    # -------------------------
    # edge 1–2 (λ3 = 0)
    for i in 1:PD - 1
        offset += 1
        L1 = _simplex_lagrange_1d(PD - i, λ1, PD)
        L2 = _simplex_lagrange_1d(i,      λ2, PD)
        dL1 = _simplex_lagrange_1d_derivative(PD - i, λ1, PD)
        dL2 = _simplex_lagrange_1d_derivative(i,      λ2, PD)

        dN[:, offset] =
            dL1 * L2 * gλ1 +
            L1 * dL2 * gλ2
    end

    # edge 2–3 (λ1 = 0)
    for i in 1:PD - 1
        offset += 1
        L2 = _simplex_lagrange_1d(PD - i, λ2, PD)
        L3 = _simplex_lagrange_1d(i,      λ3, PD)
        dL2 = _simplex_lagrange_1d_derivative(PD - i, λ2, PD)
        dL3 = _simplex_lagrange_1d_derivative(i,      λ3, PD)

        dN[:, offset] =
            dL2 * L3 * gλ2 +
            L2 * dL3 * gλ3
    end

    # edge 3–1 (λ2 = 0)
    for i in 1:PD - 1
        offset += 1
        L3 = _simplex_lagrange_1d(PD - i, λ3, PD)
        L1 = _simplex_lagrange_1d(i,      λ1, PD)
        dL3 = _simplex_lagrange_1d_derivative(PD - i, λ3, PD)
        dL1 = _simplex_lagrange_1d_derivative(i,      λ1, PD)

        dN[:, offset] =
            dL3 * L1 * gλ3 +
            L3 * dL1 * gλ1
    end

    # -------------------------
    # interior
    # -------------------------
    for i in 1:PD-2, j in 1:PD-1-i
        k = PD - i - j
        offset += 1

        L1 = _simplex_lagrange_1d(i, λ1, PD)
        L2 = _simplex_lagrange_1d(j, λ2, PD)
        L3 = _simplex_lagrange_1d(k, λ3, PD)

        dL1 = _simplex_lagrange_1d_derivative(i, λ1, PD)
        dL2 = _simplex_lagrange_1d_derivative(j, λ2, PD)
        dL3 = _simplex_lagrange_1d_derivative(k, λ3, PD)

        dN[:, offset] =
            dL1 * L2 * L3 * gλ1 +
            L1 * dL2 * L3 * gλ2 +
            L1 * L2 * dL3 * gλ3
    end

    return dN
end

function shape_function_hessian(e::Tri{Lagrange, PD}, _, ξ) where PD
    λ1 = 1.0 - ξ[1] - ξ[2]
    λ2 = ξ[1]
    λ3 = ξ[2]

    H = Array{eltype(ξ), 3}(undef, 2, 2, num_cell_dofs(e))

    # barycentric gradients
    gλ1 = [-1., -1.]
    gλ2 = [1., 0.]
    gλ3 = [0., 1.]

    # outer products
    ⊗(a,b) = a * transpose(b)

    # -------------------------
    # vertices
    # -------------------------
    H[:, :, 1] = _simplex_lagrange_1d_second_derivative(PD, λ1, PD) * (gλ1 ⊗ gλ1)
    H[:, :, 2] = _simplex_lagrange_1d_second_derivative(PD, λ2, PD) * (gλ2 ⊗ gλ2)
    H[:, :, 3] = _simplex_lagrange_1d_second_derivative(PD, λ3, PD) * (gλ3 ⊗ gλ3)

    offset = 3

    # -------------------------
    # edges
    # -------------------------
    # edge 1–2
    for i in 1:PD - 1
        offset += 1

        L1  = _simplex_lagrange_1d(PD - i, λ1, PD)
        L2  = _simplex_lagrange_1d(i,      λ2, PD)
        dL1 = _simplex_lagrange_1d_derivative(PD - i, λ1, PD)
        dL2 = _simplex_lagrange_1d_derivative(i,      λ2, PD)
        d2L1 = _simplex_lagrange_1d_second_derivative(PD - i, λ1, PD)
        d2L2 = _simplex_lagrange_1d_second_derivative(i,      λ2, PD)

        H[:, :, offset] =
            d2L1 * L2 * (gλ1 ⊗ gλ1) +
            L1 * d2L2 * (gλ2 ⊗ gλ2) +
            dL1 * dL2 * (gλ1 ⊗ gλ2 + gλ2 ⊗ gλ1)
    end

    # edge 2–3
    for i in 1:PD - 1
        offset += 1

        L2  = _simplex_lagrange_1d(PD - i, λ2, PD)
        L3  = _simplex_lagrange_1d(i,      λ3, PD)
        dL2 = _simplex_lagrange_1d_derivative(PD - i, λ2, PD)
        dL3 = _simplex_lagrange_1d_derivative(i,      λ3, PD)
        d2L2 = _simplex_lagrange_1d_second_derivative(PD - i, λ2, PD)
        d2L3 = _simplex_lagrange_1d_second_derivative(i,      λ3, PD)

        H[:, :, offset] =
            d2L2 * L3 * (gλ2 ⊗ gλ2) +
            L2 * d2L3 * (gλ3 ⊗ gλ3) +
            dL2 * dL3 * (gλ2 ⊗ gλ3 + gλ3 ⊗ gλ2)
    end

    # edge 3–1
    for i in 1:PD - 1
        offset += 1

        L3  = _simplex_lagrange_1d(PD - i, λ3, PD)
        L1  = _simplex_lagrange_1d(i,      λ1, PD)
        dL3 = _simplex_lagrange_1d_derivative(PD - i, λ3, PD)
        dL1 = _simplex_lagrange_1d_derivative(i,      λ1, PD)
        d2L3 = _simplex_lagrange_1d_second_derivative(PD - i, λ3, PD)
        d2L1 = _simplex_lagrange_1d_second_derivative(i,      λ1, PD)

        H[:, :, offset] =
            d2L3 * L1 * (gλ3 ⊗ gλ3) +
            L3 * d2L1 * (gλ1 ⊗ gλ1) +
            dL3 * dL1 * (gλ3 ⊗ gλ1 + gλ1 ⊗ gλ3)
    end

    # -------------------------
    # interior
    # -------------------------
    for i in 1:PD - 2, j in 1:PD - 1 - i
        k = PD - i - j
        offset += 1

        L1  = _simplex_lagrange_1d(i, λ1, PD)
        L2  = _simplex_lagrange_1d(j, λ2, PD)
        L3  = _simplex_lagrange_1d(k, λ3, PD)

        dL1 = _simplex_lagrange_1d_derivative(i, λ1, PD)
        dL2 = _simplex_lagrange_1d_derivative(j, λ2, PD)
        dL3 = _simplex_lagrange_1d_derivative(k, λ3, PD)

        d2L1 = _simplex_lagrange_1d_second_derivative(i, λ1, PD)
        d2L2 = _simplex_lagrange_1d_second_derivative(j, λ2, PD)
        d2L3 = _simplex_lagrange_1d_second_derivative(k, λ3, PD)

        H[:, :, offset] =
            d2L1 * L2 * L3 * (gλ1 ⊗ gλ1) +
            L1 * d2L2 * L3 * (gλ2 ⊗ gλ2) +
            L1 * L2 * d2L3 * (gλ3 ⊗ gλ3) +
            dL1 * dL2 * L3 * (gλ1 ⊗ gλ2 + gλ2 ⊗ gλ1) +
            dL1 * L2 * dL3 * (gλ1 ⊗ gλ3 + gλ3 ⊗ gλ1) +
            L1 * dL2 * dL3 * (gλ2 ⊗ gλ3 + gλ3 ⊗ gλ2)
    end

    return H
end


# function shape_function_value(e::Tri{Lagrange, PD}, X, ξ) where PD
#     # barycentric coordinates
#     λ1 = 1. - ξ[1] - ξ[2]
#     λ2 = ξ[1]
#     λ3 = ξ[2]

#     N = Vector{eltype(ξ)}(undef, num_cell_dofs(e))

#     # vertices
#     # N[1] = _barycentric(PD, λ1, PD)
#     # N[2] = _barycentric(PD, λ2, PD)
#     # N[3] = _barycentric(PD, λ3, PD)
#     N[1] = _lagrange_simplex(PD, 0, 0, λ1, λ2, λ3, PD)
#     N[2] = _lagrange_simplex(0, PD, 0, λ1, λ2, λ3, PD)
#     N[3] = _lagrange_simplex(0, 0, PD, λ1, λ2, λ3, PD)

#     # -------------------------
#     # edges
#     # -------------------------
#     offset = 3
#     # edge 1–2 (λ3 = 0)
#     for i in 1:PD-1
#         offset += 1
#         N[offset] = _lagrange_simplex(PD - i, i, 0, λ1, λ2, λ3, PD)
#     end

#     # edge 2–3 (λ1 = 0)
#     for i in 1:PD-1
#         offset += 1
#         N[offset] = _lagrange_simplex(0, PD - i, i, λ1, λ2, λ3, PD)
#     end

#     # edge 3–1 (λ2 = 0)
#     for i in 1:PD-1
#         offset += 1
#         N[offset] = _lagrange_simplex(i, 0, PD - i, λ1, λ2, λ3, PD)
#     end

#     # -------------------------
#     # interior
#     # -------------------------
#     for i in 1:PD-2, j in 1:PD-1-i
#         k = PD - i - j
#         offset += 1
#         N[offset] = _lagrange_simplex(i, j, k, λ1, λ2, λ3, PD)
#     end

#     return N
# end

