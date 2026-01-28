"""
$(TYPEDEF)
"""
struct Hex{PT, PD} <: AbstractHex{PT, PD}
end

function boundary_dofs(e::Hex{Lagrange, PD}) where PD
    # 6 faces
    nv = num_vertices_per_cell(e)  # 8
    ne = num_edges(e)              # 12
    faces = face_vertices(e)       # 4 × 6
    edges = edge_vertices(e)       # 2 × 12

    # Prepare container
    # max DOFs per face: 4 vertices + 4 * (PD - 1) edge DOFs + (PD - 1)^2 face interior
    max_rows = 4 + 4 * (PD - 1) + (PD - 1)^2
    F = zeros(Int, max_rows, num_faces(e))

    # Offset counters
    edge_offset = nv + 1          # first edge DOF index
    face_offset = nv + ne * (PD - 1) + 1  # first interior face DOF index

    if PD > 1
        for f = 1:num_faces(e)
            # corner vertices
            F[1:4, f] .= faces[:, f]

            row = 5

            # edges of this face: the edges are assumed to follow vertex order
            face_edge_pairs = [
                (faces[1, f], faces[2, f]),
                (faces[2, f], faces[3, f]),
                (faces[3, f], faces[4, f]),
                (faces[4, f], faces[1, f])
            ]

            for (v1, v2) in face_edge_pairs
                # find the edge index
                idx = findfirst(ei -> 
                    (edges[1,ei] == v1 && edges[2,ei] == v2) ||
                    (edges[1,ei] == v2 && edges[2,ei] == v1), 
                    1:ne
                )
                if PD > 1
                    F[row:row + PD - 2, f] .= edge_offset:edge_offset + PD - 2
                    edge_offset += PD - 1
                    row += PD - 1
                end
            end

            # face interior DOFs
            if PD > 1
                n_interior = (PD - 1)^2
                F[row:row + n_interior - 1, f] .= face_offset:face_offset + n_interior - 1
                face_offset += n_interior
            end
        end

        # Trim unused rows
        # return F[1:row + (PD>1 ? n_interior : 0) - 1, :]
        return F
    else
        return faces
    end
end
function dof_coordinates(e::Hex{Lagrange, PD}) where PD
    if PD == 0
        return zeros(3, 1)
    end

    Xv = vertex_coordinates(e)  # 3 × 8
    coords = copy(Xv)           # start with vertices

    if PD > 1
        # --------------------------
        # edge DOFs
        # --------------------------
        edge_pts_1d = range(-1.0, 1.0, PD + 1)[2:end - 1]  # interior points along [-1,1]

        for (v1, v2) in eachcol(edge_vertices(e))
            for ξ in edge_pts_1d
                pt = (1.0 - ξ) / 2 * Xv[:, v1] + 
                     (1.0 + ξ) / 2 * Xv[:, v2]
                coords = hcat(coords, pt)
            end
        end

        # --------------------------
        # face interior DOFs
        # --------------------------
        face_pts_1d = edge_pts_1d  # same 1D points along face
        for f in 1:num_faces(e)
            v = face_vertices(e)[:, f]
            v1, v2, v3, v4 = v
            for i in face_pts_1d, j in face_pts_1d
                # bilinear interpolation in the face
                pt = (1.0 - i) * (1.0 - j) / 4 * Xv[:, v1] +
                     (1.0 + i) * (1.0 - j) / 4 * Xv[:, v2] +
                     (1.0 + i) * (1.0 + j) / 4 * Xv[:, v3] +
                     (1.0 - i) * (1.0 + j) / 4 * Xv[:, v4]
                coords = hcat(coords, pt)
            end
        end

        # --------------------------
        # cell interior DOFs
        # --------------------------
        pts_1d = edge_pts_1d
        for i in pts_1d, j in pts_1d, k in pts_1d
            # trilinear interpolation inside the hex
            pt = zeros(eltype(coords), 3)
            for (α, β, γ, vtx) in (
                (-1, -1, -1, 1), (1, -1, -1, 2), 
                (1, 1, -1, 3), (-1, 1, -1, 4),
                (-1, -1, 1, 5), (1, -1, 1, 6), 
                (1, 1, 1, 7), (-1, 1, 1, 8)
            )
                w = (1 + α * i) / 2 * (1 + β * j) / 2 * (1 + γ * k) / 2
                pt += w * Xv[:, vtx]
            end
            coords = hcat(coords, pt)
        end
    end

    return coords
end
function interior_dofs(e::Hex{Lagrange, PD}) where PD
    if PD < 2
        return Int[]
    else
        offset = 8 + 12 * (PD - 1) + 6 * (PD - 1)^2 + 1
        return offset:offset + num_interior_dofs(e)
    end
end
num_cell_dofs(::Hex{Lagrange, PD}) where PD = (PD + 1)^3
function num_interior_dofs(::Hex{Lagrange, PD}) where PD
    if PD == 0
        return 1
    elseif PD == 1
        return 0
    else
        return (PD - 1) * (PD - 1) * (PD - 1)
    end
end

function cell_quadrature_points_and_weights(e::AbstractHex, q_rule::GaussLobattoLegendre)
    ξs, ws = cell_quadrature_points_and_weights(boundary_element(boundary_element(e, 0), 0), q_rule)
    ξ_return = Matrix{eltype(ξs)}(undef, 3, length(ξs) * length(ξs) * length(ξs) * length(ξs))
    w_return = Vector{eltype(ξs)}(undef, length(ξs) * length(ξs) * length(ξs) * length(ξs))
    for (q, ξ) in enumerate(Base.Iterators.product(ξs, ξs, ξs))
        ξ_return[1, q] = ξ[1]
        ξ_return[2, q] = ξ[2]
        ξ_return[3, q] = ξ[3]
    end
    for (q, w) in enumerate(Base.Iterators.product(ws, ws, ws))
        w_return[q] = w[1] * w[2] * w[3]
    end
    return ξ_return, w_return
end

function surface_quadrature_points_and_weights(e::AbstractHex, q_rule::GaussLobattoLegendre)
    ξs, ws = cell_quadrature_points_and_weights(boundary_element(e, 0), q_rule)
  
    ξ_return = zeros(3, length(ws), 6)
    w_return = zeros(length(ws), 6)

    ξ_return[1:2, :, 1] .= ξs
    ξ_return[3, :, 1]   .= -1.
    #
    ξ_return[1, :, 2]   .= 1.
    ξ_return[2:3, :, 2] .= ξs
    #
    ξ_return[1:2, :, 3] .= ξs
    ξ_return[3, :, 3]   .= 1.
    #
    ξ_return[1, :, 4]   .= -1.
    ξ_return[2:3, :, 4] .= ξs
    #
    ξ_return[1, :, 5]   .= ξs[1, :]
    ξ_return[2, :, 5]   .= -1.
    ξ_return[3, :, 5]   .= ξs[2, :]
    #
    ξ_return[1, :, 5]   .= ξs[1, :]
    ξ_return[2, :, 5]   .= 1.
    ξ_return[3, :, 5]   .= ξs[2, :]
    #
    #
    # ξ_return[1, :, 2] .= 1.
    # ξ_return[2, :, 2] .= ξs[1, :]
    # ξ_return[1, :, 3] .= ξs[1, :]
    # ξ_return[2, :, 3] .= 1.
    # ξ_return[1, :, 4] .= -1.
    # ξ_return[2, :, 4] .= ξs[1, :]

    for n in 1:6
        w_return[:, n] .= ws
    end
    return ξ_return, w_return
end

function _edge_dof_indices(v1::Int, v2::Int, n::Int, PD::Int)
    # Number of vertices
    nv = 8

    # Check valid edge node
    @assert 1 ≤ n ≤ PD-1

    # Edges are ordered as in the standard Hex:
    # edge 1: 1→2
    # edge 2: 2→3
    # edge 3: 3→4
    # edge 4: 4→1
    # edge 5: 5→6
    # edge 6: 6→7
    # edge 7: 7→8
    # edge 8: 8→5
    # edge 9: 1→5
    # edge10:2→6
    # edge11:3→7
    # edge12:4→8

    # List of edges in vertex ordering
    edges = [
        (1,2), (2,3), (3,4), (4,1),
        (5,6), (6,7), (7,8), (8,5),
        (1,5), (2,6), (3,7), (4,8)
    ]

    # find edge number
    edge_num = findfirst(e -> e == (v1, v2) || e == (v2, v1), edges)
    @assert edge_num !== nothing "Edge not found"

    # Each edge has PD-1 interior DOFs
    edge_offset = nv + (edge_num-1)*(PD-1)

    # global DOF index
    return edge_offset + n
end

function _face_dof_indices(face::Int, i::Int, j::Int, PD::Int)
    # i,j in 2:PD (interior nodes)
    # returns (xi_idx, eta_idx, zeta_idx)
    if face == 1          # bottom, ζ = 1
        return (i, j, 1)
    elseif face == 2      # top, ζ = PD+1
        return (i, j, PD+1)
    elseif face == 3      # left, ξ = 1
        return (1, i, j)
    elseif face == 4      # right, ξ = PD+1
        return (PD+1, i, j)
    elseif face == 5      # front, η = 1
        return (i, 1, j)
    elseif face == 6      # back, η = PD+1
        return (i, PD+1, j)
    else
        error("Invalid face index $face")
    end
end

function shape_function_value(::Hex{Lagrange, 0}, _, _)
    return ones(1)
end

function shape_function_gradient(::Hex{Lagrange, 0}, _, _)
    return zeros(1, 3)
end

function shape_function_hessian(::Hex{Lagrange, 0}, _, _)
    return zeros(1, 3, 3)
end

function shape_function_value(::Hex{Lagrange, 1}, _, ξ)
    Ns = [
        0.125 * (1 - ξ[1]) * (1 - ξ[2]) * (1 - ξ[3]),
        0.125 * (1 + ξ[1]) * (1 - ξ[2]) * (1 - ξ[3]),
        0.125 * (1 + ξ[1]) * (1 + ξ[2]) * (1 - ξ[3]),
        0.125 * (1 - ξ[1]) * (1 + ξ[2]) * (1 - ξ[3]),
        0.125 * (1 - ξ[1]) * (1 - ξ[2]) * (1 + ξ[3]),
        0.125 * (1 + ξ[1]) * (1 - ξ[2]) * (1 + ξ[3]),
        0.125 * (1 + ξ[1]) * (1 + ξ[2]) * (1 + ξ[3]),
        0.125 * (1 - ξ[1]) * (1 + ξ[2]) * (1 + ξ[3]),
    ]
    return Ns
end

function shape_function_gradient(::Hex{Lagrange, 1}, _, ξ)
    Ns = zeros(8, 3)
    # Ns = reshape([
    #   -0.125 * (1 - ξ[2]) * (1 - ξ[3]),
    #   0.125 * (1 - ξ[2]) * (1 - ξ[3]),
    #   0.125 * (1 + ξ[2]) * (1 - ξ[3]),
    #   -0.125 * (1 + ξ[2]) * (1 - ξ[3]),
    #   -0.125 * (1 - ξ[2]) * (1 + ξ[3]),
    #   0.125 * (1 - ξ[2]) * (1 + ξ[3]),
    #   0.125 * (1 + ξ[2]) * (1 + ξ[3]),
    #   -0.125 * (1 + ξ[2]) * (1 + ξ[3]),
    #   #
    #   -0.125 * (1 - ξ[1]) * (1 - ξ[3]),
    #   -0.125 * (1 + ξ[1]) * (1 - ξ[3]),
    #   0.125 * (1 + ξ[1]) * (1 - ξ[3]),
    #   0.125 * (1 - ξ[1]) * (1 - ξ[3]),
    #   -0.125 * (1 - ξ[1]) * (1 + ξ[3]),
    #   -0.125 * (1 + ξ[1]) * (1 + ξ[3]),
    #   0.125 * (1 + ξ[1]) * (1 + ξ[3]),
    #   0.125 * (1 - ξ[1]) * (1 + ξ[3]),
    #   #
    #   -0.125 * (1 - ξ[1]) * (1 - ξ[2]),
    #   -0.125 * (1 + ξ[1]) * (1 - ξ[2]),
    #   -0.125 * (1 + ξ[1]) * (1 + ξ[2]),
    #   -0.125 * (1 - ξ[1]) * (1 + ξ[2]),
    #   0.125 * (1 - ξ[1]) * (1 - ξ[2]),
    #   0.125 * (1 + ξ[1]) * (1 - ξ[2]),
    #   0.125 * (1 + ξ[1]) * (1 + ξ[2]),
    #   0.125 * (1 - ξ[1]) * (1 + ξ[2])
    # ], 8, 3)' |> collect
    Ns[1, 1] = -0.125 * (1 - ξ[2]) * (1 - ξ[3])
    Ns[2, 1] = 0.125 * (1 - ξ[2]) * (1 - ξ[3])
    Ns[3, 1] = 0.125 * (1 + ξ[2]) * (1 - ξ[3])
    Ns[4, 1] = -0.125 * (1 + ξ[2]) * (1 - ξ[3])
    Ns[5, 1] = -0.125 * (1 - ξ[2]) * (1 + ξ[3])
    Ns[6, 1] = 0.125 * (1 - ξ[2]) * (1 + ξ[3])
    Ns[7, 1] = 0.125 * (1 + ξ[2]) * (1 + ξ[3])
    Ns[8, 1] = -0.125 * (1 + ξ[2]) * (1 + ξ[3])
    #
    Ns[1, 2] = -0.125 * (1 - ξ[1]) * (1 - ξ[3])
    Ns[2, 2] = -0.125 * (1 + ξ[1]) * (1 - ξ[3])
    Ns[3, 2] = 0.125 * (1 + ξ[1]) * (1 - ξ[3])
    Ns[4, 2] = 0.125 * (1 - ξ[1]) * (1 - ξ[3])
    Ns[5, 2] = -0.125 * (1 - ξ[1]) * (1 + ξ[3])
    Ns[6, 2] = -0.125 * (1 + ξ[1]) * (1 + ξ[3])
    Ns[7, 2] = 0.125 * (1 + ξ[1]) * (1 + ξ[3])
    Ns[8, 2] = 0.125 * (1 - ξ[1]) * (1 + ξ[3])
    #
    Ns[1, 3] = -0.125 * (1 - ξ[1]) * (1 - ξ[2])
    Ns[2, 3] = -0.125 * (1 + ξ[1]) * (1 - ξ[2])
    Ns[3, 3] = -0.125 * (1 + ξ[1]) * (1 + ξ[2])
    Ns[4, 3] = -0.125 * (1 - ξ[1]) * (1 + ξ[2])
    Ns[5, 3] = 0.125 * (1 - ξ[1]) * (1 - ξ[2])
    Ns[6, 3] = 0.125 * (1 + ξ[1]) * (1 - ξ[2])
    Ns[7, 3] = 0.125 * (1 + ξ[1]) * (1 + ξ[2])
    Ns[8, 3] = 0.125 * (1 - ξ[1]) * (1 + ξ[2])
    return Ns
end

function shape_function_hessian(::Hex{Lagrange, 1}, _, _)
    return zeros(3, 3, 8)
end

function shape_function_value(e::Hex{Lagrange, PD}, _, ξ) where PD
    # ξ = [ξ, η, ζ] in reference coordinates

    le = boundary_element(boundary_element(e))
    # 1D reference shape functions along [0..PD] nodes
    coords_1d = dof_coordinates(le)  # 1D coords along [-1,1]
    N_1D = shape_function_value(le, coords_1d, ξ[1])
    M_1D = shape_function_value(le, coords_1d, ξ[2])
    L_1D = shape_function_value(le, coords_1d, ξ[3])

    N = Vector{eltype(ξ)}(undef, num_cell_dofs(e))
    idx = 1

    # ------------------------
    # 1) vertices (corners)
    # ------------------------
    for v in 1:8
        # vertex mapping: standard Hex
        # 1: (-1,-1,-1), 2:(1,-1,-1), 3:(1,1,-1), 4:(-1,1,-1)
        # 5: (-1,-1,1), etc.
        # The vertex index corresponds to 1D index along ξ, η, ζ
        i = (v == 1 || v == 4 || v == 5 || v == 8) ? 1 : PD+1
        j = (v == 1 || v == 2 || v == 5 || v == 6) ? 1 : PD+1
        k = (v <= 4) ? 1 : PD+1
        N[idx] = N_1D[i] * M_1D[j] * L_1D[k]
        idx += 1
    end

    # ------------------------
    # 2) edges
    # ------------------------
    # 12 edges: each edge varies along one coordinate with others fixed
    # Use edge_vertices(e) to know which vertices it connects
    for (v1, v2) in eachcol(edge_vertices(e))
        # find 1D index along edge direction
        for n in 2:PD
            # ξ along edge: linear mapping from v1->v2
            # determine which 1D index varies
            # The other 2 coordinates are fixed to vertex v1
            i1 = (v1 <= 4 ? 1 : PD+1)
            j1 = (v1==1 || v1==2 || v1==5 || v1==6 ? 1 : PD+1)
            k1 = (v1 <= 4 ? 1 : PD+1)
            # identify which coordinate varies along edge
            if v2 - v1 in (1,5)  # ξ varies
                xi = n
                yi = j1
                zi = k1
            elseif v2 - v1 in (2,6) # η varies
                xi = i1
                yi = n
                zi = k1
            else # ζ varies
                xi = i1
                yi = j1
                zi = n
            end
            N[idx] = N_1D[xi] * M_1D[yi] * L_1D[zi]
            idx += 1
        end
    end

    # ------------------------
    # 3) face interiors
    # ------------------------
    # 6 faces: each face varies along 2 coordinates, 1 fixed
    for f in 1:num_faces(e)
        verts = face_vertices(e)[:, f]
        # determine fixed coordinate (ξ,η,ζ)
        # use 1D indices for interior points 2:PD
        for i in 2:PD
            for j in 2:PD
                # map i,j to the two varying directions depending on face
                # fixed coordinate = first vertex's fixed direction
                # Example: bottom face ξ-η plane at ζ=1
                # This part requires careful face->coordinate mapping
                # For simplicity assume function face_dof_indices(f, i, j) exists
                xi, yi, zi = _face_dof_indices(f, i, j, PD)
                N[idx] = N_1D[xi] * M_1D[yi] * L_1D[zi]
                idx += 1
            end
        end
    end

    # ------------------------
    # 4) interior cell DOFs
    # ------------------------
    if PD > 2
        for i in 2:PD
            for j in 2:PD
                for k in 2:PD
                    N[idx] = N_1D[i] * M_1D[j] * L_1D[k]
                    idx += 1
                end
            end
        end
    end

    return N
end

function shape_function_gradient(e::Hex{Lagrange, PD}, _, ξ) where PD
    le = boundary_element(boundary_element(e))

    coords_1d = dof_coordinates(le)

    Nx  = shape_function_value(le, coords_1d, ξ[1])
    Ny  = shape_function_value(le, coords_1d, ξ[2])
    Nz  = shape_function_value(le, coords_1d, ξ[3])

    dNx = shape_function_gradient(le, coords_1d, ξ[1])
    dNy = shape_function_gradient(le, coords_1d, ξ[2])
    dNz = shape_function_gradient(le, coords_1d, ξ[3])

    ndofs = num_cell_dofs(e)
    G = Matrix{eltype(ξ)}(undef, ndofs, 3)

    idx = 1

    # -----------------------
    # vertices
    # -----------------------
    for k in (1, PD + 1), j in (1, PD + 1), i in (1, PD + 1)
        G[idx, 1] = dNx[i] * Ny[j] * Nz[k]
        G[idx, 2] = Nx[i] * dNy[j] * Nz[k]
        G[idx, 3] = Nx[i] * Ny[j] * dNz[k]
        idx += 1
    end

    # -----------------------
    # edges
    # -----------------------
    for (v1, v2) in eachcol(boundary_dofs(e))
        for n in 2:PD
            i, j, k = _edge_dof_indices(v1, v2, n, PD)

            G[idx, 1] = dNx[i] * Ny[j] * Nz[k]
            G[idx, 2] = Nx[i] * dNy[j] * Nz[k]
            G[idx, 3] = Nx[i] * Ny[j] * dNz[k]
            idx += 1
        end
    end

    # -----------------------
    # face interiors
    # -----------------------
    for f in 1:6
        for iₗ in 2:PD, jₗ in 2:PD
            i, j, k = _face_dof_indices(f, iₗ, jₗ, PD)

            G[idx, 1] = dNx[i] * Ny[j] * Nz[k]
            G[idx, 2] = Nx[i] * dNy[j] * Nz[k]
            G[idx, 3] = Nx[i] * Ny[j] * dNz[k]
            idx += 1
        end
    end

    # -----------------------
    # cell interior
    # -----------------------
    if PD > 2
        for i in 2:PD, j in 2:PD, k in 2:PD
            G[idx, 1] = dNx[i] * Ny[j] * Nz[k]
            G[idx, 2] = Nx[i] * dNy[j] * Nz[k]
            G[idx, 3] = Nx[i] * Ny[j] * dNz[k]
            idx += 1
        end
    end

    return G
end

# function shape_function_hessian(e::Hex{Lagrange, PD}, _, ξ) where PD
#     le = boundary_element(boundary_element(e))

#     coords_1d = dof_coordinates(le)
#     Nx   = shape_function_value(le, coords_1d, ξ[1])
#     Ny   = shape_function_value(le, coords_1d, ξ[2])
#     Nz   = shape_function_value(le, coords_1d, ξ[3])

#     dNx  = shape_function_gradient(le, coords_1d, ξ[1])
#     dNy  = shape_function_gradient(le, coords_1d, ξ[2])
#     dNz  = shape_function_gradient(le, coords_1d, ξ[3])

#     d2Nx = shape_function_hessian(le, coords_1d, ξ[1])
#     d2Ny = shape_function_hessian(le, coords_1d, ξ[2])
#     d2Nz = shape_function_hessian(le, coords_1d, ξ[3])

#     ndofs = num_cell_dofs(e)
#     H = Array{eltype(ξ)}(undef, 3, 3, ndofs)

#     idx = 1

#     # -----------------------
#     # vertices
#     # -----------------------
#     for k in (1, PD + 1), j in (1, PD + 1), i in (1, PD + 1)
#         H[:, :, idx] = [
#             d2Nx[i] * Ny[j] * Nz[k] dNx[i] * dNy[j] * Nz[k] dNx[i] * Ny[j] * dNz[k];
#             dNx[i] * dNy[j] * Nz[k] Nx[i] * d2Ny[j] * Nz[k] Nx[i] * dNy[j] * dNz[k];
#             dNx[i] * Ny[j] * dNz[k] Nx[i] * dNy[j] * dNz[k] Nx[i] * Ny[j] * d2Nz[k]
#         ]
#         idx += 1
#     end

#     # -----------------------
#     # edges
#     # -----------------------
#     for (v1, v2) in eachcol(edge_vertices(e))
#         for n in 2:PD
#             i, j, k = edge_dof_indices(v1, v2, n, PD)
#             H[:, :, idx] = [
#                 d2Nx[i] * Ny[j] * Nz[k] dNx[i] * dNy[j] * Nz[k] dNx[i] * Ny[j] * dNz[k];
#                 dNx[i] * dNy[j] * Nz[k] Nx[i] * d2Ny[j] * Nz[k] Nx[i] * dNy[j] * dNz[k];
#                 dNx[i] * Ny[j] * dNz[k] Nx[i] * dNy[j] * dNz[k] Nx[i] * Ny[j] * d2Nz[k]
#             ]
#             idx += 1
#         end
#     end

#     # -----------------------
#     # face interiors
#     # -----------------------
#     for f in 1:6
#         for il in 2:PD, jl in 2:PD
#             i, j, k = face_dof_indices(f, il, jl, PD)

#             H[:, :, idx] = [
#                 d2Nx[i] * Ny[j] * Nz[k] dNx[i] * dNy[j] * Nz[k] dNx[i] * Ny[j] * dNz[k];
#                 dNx[i] * dNy[j] * Nz[k] Nx[i] * d2Ny[j] * Nz[k] Nx[i] * dNy[j] * dNz[k];
#                 dNx[i] * Ny[j] * dNz[k] Nx[i] * dNy[j] * dNz[k] Nx[i] * Ny[j] * d2Nz[k]
#             ]
#             idx += 1
#         end
#     end

#     # -----------------------
#     # cell interior
#     # -----------------------
#     if PD > 2
#         for i in 2:PD, j in 2:PD, k in 2:PD
#             H[:, :, idx] = [
#                 d2Nx[i] * Ny[j] * Nz[k] dNx[i] * dNy[j] * Nz[k] dNx[i] * Ny[j] * dNz[k];
#                 dNx[i] * dNy[j] * Nz[k] Nx[i] * d2Ny[j] * Nz[k] Nx[i] * dNy[j] * dNz[k];
#                 dNx[i] * Ny[j] * dNz[k] Nx[i] * dNy[j] * dNz[k] Nx[i] * Ny[j] * d2Nz[k]
#             ]
#             idx += 1
#         end
#     end

#     return H
# end
