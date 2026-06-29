"""
$(TYPEDEF)
"""
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

function interior_dofs(::Tri{Lagrange, PD}) where PD
    if PD < 3
        return Int[]
    else
        @assert false "TODO"
    end
end

num_cell_dofs(::Tri{Lagrange, PD}) where PD = (PD + 1) * (PD + 2) ÷ 2
num_dofs_on_boundary(::Tri{Lagrange, PD}, ::Int) where PD = PD == 0 ? 2 : PD + 1
num_interior_dofs(::Tri{Lagrange, PD}) where PD = PD < 3 ? 0 : (PD - 1) * (PD - 2) ÷ 2

function shape_function_value(::Tri{Lagrange, 0}, _)
    return ones(1)
end

function shape_function_gradient(::Tri{Lagrange, 0}, _)
    return zeros(1, 2)
end

function shape_function_hessian(::Tri{Lagrange, 0}, _)
    return zeros(1, 2, 2)
end

function shape_function_value(::Tri{Lagrange, PD}, ξ) where PD
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

function shape_function_gradient(e::Tri{Lagrange, PD}, ξ) where PD
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

function shape_function_hessian(e::Tri{Lagrange, PD}, ξ) where PD
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

function geometry_shape_function_value(::Tri{RaviartThomas, 0}, ξ)
    return shape_function_value(Tri{Lagrange, 1}(), ξ)
end

function geometry_shape_function_gradient(::Tri{RaviartThomas, 0}, ξ)
    return shape_function_gradient(Tri{Lagrange, 1}(), ξ)
end

# https://defelement.org/elements/examples/triangle-raviart-thomas-lagrange-0.html
# but re-ordered for exodus numbering
function shape_function_value(::Tri{RaviartThomas, 0}, ξ)
    N = Matrix{Float64}(undef, 3, 2)
    #
    N[1, 1] = -ξ[1]
    N[1, 2] = 1 - ξ[2]
    #
    N[2, 1] = -ξ[1]
    N[2, 2] = -ξ[2]
    #
    N[3, 1] = ξ[1] - 1
    N[3, 2] = ξ[2]

    return N
end

function shape_function_divergence(::Tri{RaviartThomas, 0}, ξ)
    return [-2., -2., 2.]
end
