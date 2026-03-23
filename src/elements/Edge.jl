const EQUALLY_SPACED = 1
const GLL            = 2

"""
$(TYPEDEF)
"""
struct Edge{PT, PD} <: AbstractEdge{PT, PD}
    nodal_locations::Int
    shifted::Bool

    function Edge{PT, PD}(; 
        nodal_locations::Int = EQUALLY_SPACED,
        shifted::Bool = false
    ) where {PT, PD}
        @assert nodal_locations == EQUALLY_SPACED ||
                nodal_locations == GLL
        return new{PT, PD}(nodal_locations, shifted)
    end

    function Edge{Hermite, PD}(; 
        nodal_locations::Int = EQUALLY_SPACED,
        shifted::Bool = false
    ) where PD
        @assert PD >= 3 "PD should be greater than or equal to 3 for Hermite"
        @assert nodal_locations == EQUALLY_SPACED ||
                nodal_locations == GLL
        return new{Hermite, PD}(nodal_locations, shifted)
    end
end

# dof methods
boundary_dofs(::Edge) = [1 2]' |> collect

# Hermite specific
# TODO finish me
function dof_coordinates(e::Edge{Hermite, PD}) where PD
    if e.shifted
        @assert false
    end
    modal = [basis(Legendre, k) for k in 0:PD]

    # DOF locations
    interior = PD > 3 ? collect(range(-1, 1; length=PD - 1)[2:end - 1]) : Float64[]
    ndofs = 4 + length(interior)

    @assert ndofs == PD + 1

    # Assemble interpolation matrix
    A = zeros(Float64, ndofs, PD + 1)

    row = 1

    # u(-1)
    for j in 1:PD + 1
        A[row, j] = modal[j](-1)
    end
    row += 1

    # u'(-1)
    for j in 1:PD + 1
        A[row, j] = derivative(modal[j])(-1)
    end
    row += 1

    # u(1)
    for j in 1:PD + 1
        A[row, j] = modal[j](1)
    end
    row += 1

    # u'(1)
    for j in 1:PD + 1
        A[row, j] = derivative(modal[j])(1)
    end
    row += 1

    # interior value DOFs
    for ξ in interior
        for j in 1:PD + 1
            A[row, j] = modal[j](ξ)
        end
        row += 1
    end

    # Solve for nodal basis
    coeffs = A \ I(ndofs)

    return [
        sum(coeffs[j, i] * modal[j] for j in 1:PD + 1)
        for i in 1:ndofs
    ]
end
interior_dofs(e::Edge{Hermite, PD}) where PD = 5:num_cell_dofs(e) |> collect
num_cell_dofs(e::Edge{Hermite, PD}) where PD = 4 + num_interior_dofs(e)
num_interior_dofs(::Edge{Hermite, PD}) where PD = PD > 3 ? PD - 1 : 0

# Lagrange specific
function dof_coordinates(e::Edge{Lagrange, PD}) where PD
    if e.shifted
        x_mid = 0.5
        x_min = 0.
    else
        x_mid = 0.
        x_min = -1.
    end

    if PD == 0
        coords = zeros(1, 1)
    elseif PD == 1
        coords = zeros(1, PD + 1)
        coords[1, 1] = x_min
        coords[1, 2] = 1.
    elseif PD == 2
        coords = zeros(1, PD + 1)
        coords[1, 1] = x_min
        coords[1, 2] = 1.
        coords[1, 3] = x_mid
    else
        coords_rng = LinRange(x_min, 1., PD + 1)
        coords = zeros(1, PD + 1)
        coords[1, 1] = coords_rng[1]
        coords[1, 2] = coords_rng[end]

        if e.nodal_locations == EQUALLY_SPACED
            coords_rng = LinRange(x_mid, 1., PD + 1)
            coords[1, 3:end] .= coords_rng[2:end - 1]
        elseif e.nodal_locations == GLL
            @assert false "Need to finish this up"
            # poly_coeffs = zeros(Int, PD + 1)
            # poly_coeffs[end] = 1
            # poly = Legendre(poly_coeffs)
            # rts = roots(poly)
            # @assert all(imag.(rts) .== 0.)
            # coords[1, 3:end] .= real(rts)
        end
    end
    return coords
end
interior_dofs(e::Edge{Lagrange, PD}) where PD = 3:num_cell_dofs(e) |> collect
num_cell_dofs(::Edge{Lagrange, PD}) where PD = PD + 1
num_dofs_on_boundary(::Edge{Lagrange, PD}, ::Int) where PD = 1
num_interior_dofs(::Edge{Lagrange, PD}) where PD = PD < 2 ? 0 : PD - 1

function cell_quadrature_points_and_weights(e::AbstractEdge, q_rule::GaussLegendre)
    ξs, ws = gausslegendre(cell_quadrature_degree(q_rule))

    if e.shifted
        ξs .= (ξs .+ 1.) ./ 2.
        ws .= ws ./ 2.
    end

    return reshape(ξs, 1, length(ξs)), ws
end

function surface_quadrature_points_and_weights(e::AbstractEdge, ::GaussLegendre)
    if e.shifted
        x_min = 0.
    else
        x_min = -1.
    end

    ξs = zeros(1, 1, 2)
    ξs[1, 1, 1] = x_min
    ξs[1, 1, 2] = 1.
    ws = ones(1, 2)
    return ξs, ws
end

function cell_quadrature_points_and_weights(e::AbstractEdge, q_rule::GaussLobattoLegendre)
    ξs, ws = gausslegendre(cell_quadrature_degree(q_rule))

    if e.shifted
        ξs .= (ξs .+ 1.) ./ 2.
        ws .= ws ./ 2.
    end

    return reshape(ξs, 1, length(ξs)), ws
end

function surface_quadrature_points_and_weights(e::AbstractEdge, ::GaussLobattoLegendre)
    if e.shifted
        x_min = 0.
    else
        x_min = -1.
    end

    ξs = zeros(1, 1, 2)
    ξs[1, 1, 1] = x_min
    ξs[1, 1, 2] = 1.
    ws = ones(1, 2)
    return ξs, ws
end

# 0th order Lagrange implementation
function shape_function_value(::Edge{Lagrange, 0}, _, ::Number)
    return ones(1)
end

function shape_function_gradient(::Edge{Lagrange, 0}, _, ::Number)
    return zeros(1, 1)
end

function shape_function_hessian(::Edge{Lagrange, 0}, _, ::Number)
    return zeros(1, 1, 1)
end

# nth order Lagrange implementation
function shape_function_value(e::Edge{Lagrange, PD}, Xs, ξ::Number) where PD
    if e.shifted
        type = ShiftedLegendre
    else
        type = Legendre
    end

    n_nodes = polynomial_degree(e) + 1
    A = zeros(length(Xs), n_nodes)
    nf = zeros(1, n_nodes)
    for n in 1:n_nodes
        p = basis(type, n - 1)
        p = sqrt(2 * (n - 1) + 1) * p
        for m in axes(A, 1)
            A[m, n] = p(Xs[m][1])
        end
        nf[1, n] = p(ξ[1])
    end
    N = A' \ nf'
    return N[:, 1]
end

function shape_function_gradient(e::Edge{Lagrange, PD}, Xs, ξ) where PD
    if e.shifted
        type = ShiftedLegendre
    else
        type = Legendre
    end

    n_nodes = polynomial_degree(e) + 1
    A = zeros(length(Xs), n_nodes)
    nf = zeros(1, n_nodes)
    for n in 1:n_nodes
        p = basis(type, n - 1)
        p = sqrt(2 * (n - 1) + 1) * p
        for m in axes(A, 1)
            A[m, n] = p(Xs[m][1])
        end
        nf[1, n] = derivative(p)(ξ[1])
    end
    ∇N_ξ = A' \ nf'
    return reshape(∇N_ξ, n_nodes, 1)
end

function shape_function_hessian(e::Edge{Lagrange, PD}, Xs, ξ) where PD
    if e.shifted
        type = ShiftedLegendre
    else
        type = Legendre
    end

    n_nodes = polynomial_degree(e) + 1
    A = zeros(length(Xs), n_nodes)
    nf = zeros(1, n_nodes)
    for n in 1:n_nodes
        p = basis(type, n - 1)
        p = sqrt(2 * (n - 1) + 1) * p
        for m in axes(A, 1)
            A[m, n] = p(Xs[m][1])
        end

        # need to handle this case specially
        # in SpecialPolynomials v4
        if n == 1
            nf[1, n] = 0.
        else
            nf[1, n] = derivative(p, 2)(ξ[1])
        end
    end
    ∇∇N_ξ = A' \ nf'
    return reshape(∇∇N_ξ, n_nodes, 1, 1)
end

# nth order Serendipity implementation
