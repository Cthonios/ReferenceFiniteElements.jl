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
function dof_coordinates(e::Quad{Lagrange, PD}) where PD
    # if PD == 0
    #     return zeros(2, 1)
    # end

    coords = vertex_coordinates(e)[1:2, :]

    if PD > 1
        # do edge midpoints
        edge_coords = dof_coordinates(boundary_element(e, 0))

        # face 1
        for n in 1:PD - 1
            coords = hcat(coords, [edge_coords[1, n + 2], -1.])
        end
        # face 2
        for n in 1:PD - 1
            coords = hcat(coords, [1., edge_coords[1, n + 2]])
        end
        # face 3
        for n in 1:PD - 1
            coords = hcat(coords, [edge_coords[1, n + 2], 1.])
        end
        # face 4
        for n in 1:PD - 1
            coords = hcat(coords, [-1., edge_coords[1, n + 2]])
        end

        # now for interiors
        for n in 1:PD - 1
            for m in 1:PD - 1
                coords = hcat(coords, [edge_coords[1, m + 2], edge_coords[1, n + 2]])
            end
        end
    end
    return coords
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
function num_interior_dofs(::Quad{Lagrange, PD}) where PD
    if PD == 0
        return 1
    elseif PD == 1
        return 0
    else
        return (PD - 1) * (PD - 1)
    end
end

function cell_quadrature_points_and_weights(e::AbstractQuad, q_rule::GaussLobattoLegendre)
    # ξs, ws = gausslegendre(cell_quadrature_degree(e))
    ξs, ws = cell_quadrature_points_and_weights(boundary_element(e, 0), q_rule)
    ξ_return = Matrix{eltype(ξs)}(undef, 2, length(ws) * length(ws))
    w_return = Vector{eltype(ξs)}(undef, length(ws) * length(ws))
    for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs))
        ξ_return[1, q] = ξ[1]
        ξ_return[2, q] = ξ[2]
    end
    for (q, w) in enumerate(Base.Iterators.product(ws, ws))
        w_return[q] = w[1] * w[2]
    end
    return ξ_return, w_return
end

function surface_quadrature_points_and_weights(e::AbstractQuad, q_rule::GaussLobattoLegendre)
    ξs, ws = cell_quadrature_points_and_weights(boundary_element(e, 0), q_rule)
  
    ξ_return = zeros(2, length(ws), 4)
    w_return = zeros(length(ws), 4)

    ξ_return[1, :, 1] .= ξs[1, :]
    ξ_return[2, :, 1] .= -1.
    ξ_return[1, :, 2] .= 1.
    ξ_return[2, :, 2] .= ξs[1, :]
    ξ_return[1, :, 3] .= ξs[1, :]
    ξ_return[2, :, 3] .= 1.
    ξ_return[1, :, 4] .= -1.
    ξ_return[2, :, 4] .= ξs[1, :]

    for n in 1:4
        w_return[:, n] .= ws
    end
    return ξ_return, w_return
end

function shape_function_value(::Quad{Lagrange, 0}, _, _)
    return ones(1)
end

function shape_function_gradient(::Quad{Lagrange, 0}, _, _)
    return zeros(1, 2)
end

function shape_function_hessian(::Quad{Lagrange, 0}, _, _)
    return zeros(1, 2, 2)
end

function shape_function_value(e::Quad{Lagrange, PD}, _, ξ) where PD
    coords_x = dof_coordinates(boundary_element(e, 0))
    coords_y = dof_coordinates(boundary_element(e, 0))
    N_x = shape_function_value(boundary_element(e, 0), coords_x, ξ[1])
    N_y = shape_function_value(boundary_element(e, 0), coords_y, ξ[2])
  
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

function shape_function_gradient(e::Quad{Lagrange, PD}, X, ξ) where PD
    coords_x = dof_coordinates(boundary_element(e, 0))
    coords_y = dof_coordinates(boundary_element(e, 0))
    N_x = shape_function_value(boundary_element(e, 0), coords_x, ξ[1])
    N_y = shape_function_value(boundary_element(e, 0), coords_y, ξ[2])
    ∇N_x = shape_function_gradient(boundary_element(e, 0), coords_x, ξ[1])
    ∇N_y = shape_function_gradient(boundary_element(e, 0), coords_y, ξ[2])
  
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

function shape_function_hessian(e::Quad{Lagrange, PD}, X, ξ) where PD
    coords_x = dof_coordinates(boundary_element(e, 0))
    coords_y = dof_coordinates(boundary_element(e, 0))
    N_x = shape_function_value(boundary_element(e, 0), coords_x, ξ[1])
    N_y = shape_function_value(boundary_element(e, 0), coords_y, ξ[2])
    ∇N_x = shape_function_gradient(boundary_element(e, 0), coords_x, ξ[1])
    ∇N_y = shape_function_gradient(boundary_element(e, 0), coords_y, ξ[2])
    ∇∇N_x = shape_function_hessian(boundary_element(e, 0), coords_x, ξ[1])
    ∇∇N_y = shape_function_hessian(boundary_element(e, 0), coords_y, ξ[2])
  
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
