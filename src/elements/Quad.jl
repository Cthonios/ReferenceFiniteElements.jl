"""
$(TYPEDEF)
"""
struct Quad{PT, PD} <: AbstractQuad{PT, PD}
end

function boundary_dofs(e::Quad{Lagrange, PD}) where PD
    linear_edges = edge_vertices(e)
    if PD < 2
        return linear_edges
    else
        edges = zeros(Int, PD + 1, 4)
        edges[1:2, 1:4] .= linear_edges
        offset = 5
        for n in 1:4
            edges[3:end, n] = offset:offset + PD - 2
            offset += PD - 1
        end
        return edges
    end
end

function interior_dofs(e::Quad{Lagrange, PD}) where PD
    if PD < 2
        return Int[]
    else
        offset = 4 + 4 * (PD - 1) + 1
        return offset:offset + num_interior_dofs(e) - 1 |> collect
    end
end

num_cell_dofs(::Quad{Lagrange, PD}) where PD = (PD + 1) * (PD + 1)
# provides bdofs...?
num_dofs_on_boundary(::Quad{Lagrange, PD}, ::Int) where PD = PD == 0 ? 2 : PD + 1

function num_interior_dofs(::Quad{Lagrange, PD}) where PD
    if PD == 0
        return 1
    elseif PD == 1
        return 0
    else
        return (PD - 1) * (PD - 1)
    end
end

function shape_function_value(::Quad{Lagrange, 0}, _)
    return ones(1)
end

function shape_function_gradient(::Quad{Lagrange, 0}, _)
    return zeros(1, 2)
end

function shape_function_hessian(::Quad{Lagrange, 0}, _)
    return zeros(1, 2, 2)
end

function shape_function_value(e::Quad{Lagrange, PD}, ξ) where PD
    N_x = shape_function_value(boundary_element(e, 0), ξ[1])
    N_y = shape_function_value(boundary_element(e, 0), ξ[2])
  
    N = Vector{eltype(ξ)}(undef, num_cell_dofs(e))
  
    # corner nodes first
    N[1] = N_x[1] * N_y[1]
    N[2] = N_x[2] * N_y[1]
    N[3] = N_x[2] * N_y[2]
    N[4] = N_x[1] * N_y[2]

    # facet 1
    offset = 4
    for n in 2:PD
        N[offset + n - 1] = N_x[n + 1] * N_y[1]
    end

    # facet 2
    offset += PD - 1
    for n in 2:PD
        N[offset + n - 1] = N_x[2] * N_y[n + 1]
    end

    # facet 3
    offset += PD - 1
    for n in 2:PD
        N[offset + n - 1] = N_x[n + 1] * N_y[2]
    end

    # facet 4
    offset += PD - 1
    for n in 2:PD
        N[offset + n - 1] = N_x[1] * N_y[n + 1]
    end
    
    # now for interior nodes
    m = 4 + 4 * (PD - 1) + 1
    for (N_1, N_2) in Iterators.product(N_x[3:end], N_y[3:end])
        N[m] = N_1 * N_2
        m = m + 1
    end 
    return N
end

function shape_function_gradient(e::Quad{Lagrange, PD}, ξ) where PD
    N_x = shape_function_value(boundary_element(e, 0), ξ[1])
    N_y = shape_function_value(boundary_element(e, 0), ξ[2])
    ∇N_x = shape_function_gradient(boundary_element(e, 0), ξ[1])
    ∇N_y = shape_function_gradient(boundary_element(e, 0), ξ[2])
  
    # return N_x * N_y
  
    ∇N = Matrix{eltype(ξ)}(undef, num_cell_dofs(e), 2)
  
    # corner nodes first
    ∇N[1, 1] = ∇N_x[1] * N_y[1]
    ∇N[1, 2] = N_x[1] * ∇N_y[1]
    ∇N[2, 1] = ∇N_x[2] * N_y[1]
    ∇N[2, 2] = N_x[2] * ∇N_y[1]
    ∇N[3, 1] = ∇N_x[2] * N_y[2]
    ∇N[3, 2] = N_x[2] * ∇N_y[2]
    ∇N[4, 1] = ∇N_x[1] * N_y[2]
    ∇N[4, 2] = N_x[1] * ∇N_y[2]
    
    # edge nodes next
    # facet 1
    offset = 4
    for n in 2:PD
        ∇N[offset + n - 1, 1] = ∇N_x[n + 1] * N_y[1]
        ∇N[offset + n - 1, 2] = N_x[n + 1] * ∇N_y[1]
    end

    # facet 2
    offset += PD - 1
    for n in 2:PD
        ∇N[offset + n - 1, 1] = ∇N_x[2] * N_y[n + 1]
        ∇N[offset + n - 1, 2] = N_x[2] * ∇N_y[n + 1]
    end

    # facet 3
    offset += PD - 1
    for n in 2:PD
        ∇N[offset + n - 1, 1] = ∇N_x[n + 1] * N_y[2]
        ∇N[offset + n - 1, 2] = N_x[n + 1] * ∇N_y[2]
    end

    # facet 4
    offset += PD - 1
    for n in 2:PD
        ∇N[offset + n - 1, 1] = ∇N_x[1] * N_y[n + 1]
        ∇N[offset + n - 1, 2] = N_x[1] * ∇N_y[n + 1]
    end
  
    # now for interior nodes
    m = 4 + 4 * (PD - 1) + 1
    Ns = Iterators.product(N_x[3:end], N_y[3:end])
    ∇Ns = Iterators.product(∇N_x[3:end], ∇N_y[3:end])
    for ((N_1, N_2), (∇N_1, ∇N_2)) in zip(Ns, ∇Ns)
        ∇N[m, 1] = ∇N_1 * N_2
        ∇N[m, 2] = N_1 * ∇N_2
        m = m + 1
    end 
  
    return ∇N
end

function shape_function_hessian(e::Quad{Lagrange, PD}, ξ) where PD
    N_x = shape_function_value(boundary_element(e, 0), ξ[1])
    N_y = shape_function_value(boundary_element(e, 0), ξ[2])
    ∇N_x = shape_function_gradient(boundary_element(e, 0), ξ[1])
    ∇N_y = shape_function_gradient(boundary_element(e, 0), ξ[2])
    ∇∇N_x = shape_function_hessian(boundary_element(e, 0), ξ[1])
    ∇∇N_y = shape_function_hessian(boundary_element(e, 0), ξ[2])
  
    ∇∇N = Array{eltype(ξ), 3}(undef, num_cell_dofs(e), 2, 2)
  
    # corner nodes first
    ∇∇N[1, 1, 1] = ∇∇N_x[1] * N_y[1]
    ∇∇N[1, 1, 2] = ∇N_x[1] * ∇N_y[1]
    ∇∇N[1, 2, 1] = ∇N_x[1] * ∇N_y[1]
    ∇∇N[1, 2, 2] = N_x[1] * ∇∇N_y[1]
    #
    ∇∇N[2, 1, 1] = ∇∇N_x[2] * N_y[1]
    ∇∇N[2, 1, 2] = ∇N_x[2] * ∇N_y[1]
    ∇∇N[2, 2, 1] = ∇N_x[2] * ∇N_y[1]
    ∇∇N[2, 2, 2] = N_x[2] * ∇∇N_y[1]
    #
    ∇∇N[3, 1, 1] = ∇∇N_x[2] * N_y[2]
    ∇∇N[3, 1, 2] = ∇N_x[2] * ∇N_y[2]
    ∇∇N[3, 2, 1] = ∇N_x[2] * ∇N_y[2]
    ∇∇N[3, 2, 2] = N_x[2] * ∇∇N_y[2]
    #
    ∇∇N[4, 1, 1] = ∇∇N_x[1] * N_y[2]
    ∇∇N[4, 1, 2] = ∇N_x[1] * ∇N_y[2]
    ∇∇N[4, 2, 1] = ∇N_x[1] * ∇N_y[2]
    ∇∇N[4, 2, 2] = N_x[1] * ∇∇N_y[2]
    
    # TODO need to fix midpoints and interiors
    # edge nodes next
    # facet 1
    offset = 4
    for n in 2:PD
        ∇∇N[offset + n - 1, 1, 1] = ∇∇N_x[n + 1] * N_y[1]
        ∇∇N[offset + n - 1, 1, 2] = ∇N_x[n + 1] * ∇N_y[1]
        ∇∇N[offset + n - 1, 2, 1] = ∇N_x[n + 1] * ∇N_y[1]
        ∇∇N[offset + n - 1, 2, 2] = N_x[n + 1] * ∇∇N_y[1]
    end

    # facet 2
    offset += PD - 1
    for n in 2:PD
        ∇∇N[offset + n - 1, 1, 1] = ∇∇N_x[2] * N_y[n + 1]
        ∇∇N[offset + n - 1, 1, 2] = ∇N_x[2] * ∇N_y[n + 1]
        ∇∇N[offset + n - 1, 2, 1] = ∇N_x[2] * ∇N_y[n + 1]
        ∇∇N[offset + n - 1, 2, 2] = N_x[2] * ∇∇N_y[n + 1]
    end

    # facet 3
    offset += PD - 1
    for n in 2:PD
        ∇∇N[offset + n - 1, 1, 1] = ∇∇N_x[n + 1] * N_y[2]
        ∇∇N[offset + n - 1, 1, 2] = ∇N_x[n + 1] * ∇N_y[2]
        ∇∇N[offset + n - 1, 2, 1] = ∇N_x[n + 1] * ∇N_y[2]
        ∇∇N[offset + n - 1, 2, 2] = N_x[n + 1] * ∇∇N_y[2]
    end

    # facet 4
    offset += PD - 1
    for n in 2:PD
        ∇∇N[offset + n - 1, 1, 1] = ∇∇N_x[1] * N_y[n + 1]
        ∇∇N[offset + n - 1, 1, 2] = ∇N_x[1] * ∇N_y[n + 1]
        ∇∇N[offset + n - 1, 2, 1] = ∇N_x[1] * ∇N_y[n + 1]
        ∇∇N[offset + n - 1, 2, 2] = N_x[1] * ∇∇N_y[n + 1]
    end
  
    # # now for interior nodes
    m = 4 + 4 * (PD - 1) + 1
    Ns = Iterators.product(N_x[3:end], N_y[3:end])
    ∇Ns = Iterators.product(∇N_x[3:end], ∇N_y[3:end])
    ∇∇Ns = Iterators.product(∇∇N_x[3:end], ∇∇N_y[3:end])

    for ((N_1, N_2), (∇N_1, ∇N_2), (∇∇N_1, ∇∇N_2)) in zip(Ns, ∇Ns, ∇∇Ns)
        ∇∇N[m, 1, 1] = ∇∇N_1 * N_2
        ∇∇N[m, 1, 2] = ∇N_1 * ∇N_2
        ∇∇N[m, 2, 1] = ∇N_1 * ∇N_2
        ∇∇N[m, 2, 2] = N_1 * ∇∇N_2
        m = m + 1
    end
  
    return ∇∇N
end

########################################################################
# Raviart-Thomas implementation
########################################################################
function boundary_dofs(::Quad{RaviartThomas, 0})
    return reshape(collect(1:4), 1, 4)
end
interior_dofs(::Quad{RaviartThomas, 0}) = Int[]
num_cell_dofs(::Quad{RaviartThomas, 0}) = 4
num_dofs_on_boundary(::Quad{RaviartThomas, 0}, ::Int) = 1
num_interior_dofs(::Quad{RaviartThomas, 0}) = 0

function geometry_shape_function_value(::Quad{RaviartThomas, 0}, ξ)
    return shape_function_value(Quad{Lagrange, 1}(), X, ξ)
end

function geometry_shape_function_gradient(::Quad{RaviartThomas, 0}, ξ)
    return shape_function_gradient(Quad{Lagrange, 1}(), ξ)
end

# https://defelement.org/elements/examples/quadrilateral-raviart-thomas-lagrange-0.html
# but re-ordered for exodus
function shape_function_value(::Quad{RaviartThomas, 0}, ξ)
    N = Matrix{Float64}(undef, 4, 2)

    # bottom
    N[1, 1] = 0.0
    N[1, 2] = (1.0 - ξ[2]) / 2.0

    # right
    N[2, 1] = -(ξ[1] + 1.0) / 2.0
    N[2, 2] = 0.0

    # top
    N[3, 1] = 0.0
    N[3, 2] = (ξ[2] + 1.0) / 2.0

    # left
    N[4, 1] = (ξ[1] - 1.0) / 2.0
    N[4, 2] = 0.0

    return N
end

function shape_function_divergence(::Quad{RaviartThomas, 0}, ξ)
    return [-0.5, -0.5, 0.5, 0.5]
end
