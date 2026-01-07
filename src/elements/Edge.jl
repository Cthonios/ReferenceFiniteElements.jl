struct Edge{PT, PD} <: AbstractEdge{PT, PD}
    shifted::Bool

    function Edge{PT, PD}(; shifted::Bool = false) where {PT, PD}
        return new{PT, PD}(shifted)
    end
end

# dof methods
boundary_dofs(::Edge{Lagrange, PD}) where PD = [1 2]' |> collect

function dof_coordinates(e::Edge{Lagrange, PD}) where PD
    if e.shifted
        x_min = 0.
    else
        x_min = -1.
    end

    if PD == 0
        coords = zeros(1, 1)
    elseif PD == 1
        coords = zeros(1, PD + 1)
        coords[1, 1] = x_min
        coords[1, 2] = 1.
    else
        coords_rng = LinRange(x_min, 1., PD + 1)
        coords = zeros(1, PD + 1)
        coords[1, 1] = coords_rng[1]
        coords[1, 2] = coords_rng[end]
        coords[1, 3:end] .= coords_rng[2:end - 1]
    end
    return coords
end

interior_dofs(e::Edge{Lagrange, PD}) where PD = 3:num_cell_dofs(e) |> collect

num_cell_dofs(::Edge{Lagrange, PD}) where PD = PD + 1
num_interior_dofs(::Edge{Lagrange, PD}) where PD = PD < 2 ? 0 : PD - 1

# surface_dof_coordinates(::Edge{Lagrange, PD}) where PD = [-1. 1.]' |> collect

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

# nth order implementation
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
    return reshape(∇N_ξ, 1, n_nodes)
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
    return reshape(∇∇N_ξ, 1, 1, n_nodes)
end
