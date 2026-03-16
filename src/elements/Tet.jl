"""
$(TYPEDEF)
"""
struct Tet{PT, PD} <: AbstractTet{PT, PD}
end

function boundary_dofs(e::Tet{Lagrange, PD}) where PD
    # Canonical linear vertices
    linear_vertices = cell_vertices(e)
    
    # Edge and face connectivity (1-based)
    # Provided by your conventions
    edge_verts = edge_vertices(e)
    face_verts = face_vertices(e)

    # Total number of faces
    nfaces = size(face_verts, 2)

    # Number of DOFs per face
    # For linear: 3 vertex DOFs
    # For PD >= 2, we include edge interiors
    # faces = zeros(Int, PD + 1 + (PD - 1) * (PD - 2) ÷ 2, nfaces)
    if PD < 2
        # faces = zeros(Int, 3, nfaces)
        return face_verts
    elseif PD == 2
        # faces = zeros(Int, 6, nfaces)
        return [
            2   1   1   1;
            3   4   2   3;
            4   3   4   2;
            6   8   5   7;
           10  10   9   6;
            9   7   8   5;
        ]
    else
        faces = zeros(Int, (PD - 1) * (PD - 2) ÷ 2, nfaces)
    end

    # --- Step 1: vertices on faces ---
    for f = 1:nfaces
        faces[1:3, f] .= face_verts[:,f]  # linear vertices of face
    end

    offset = 4  # first edge/face DOF index

    # --- Step 2: edge interior DOFs ---
    # Each edge has (PD - 1) interior DOFs
    if PD ≥ 2
        nedges = size(edge_verts,2)
        edge_dof_start = 5  # DOF index after vertices
        edge_offsets = zeros(Int, nedges)
        for eidx = 1:nedges
            edge_offsets[eidx] = offset
            offset += PD-1
        end

        # Map edge interiors to faces
        for f = 1:nfaces
            for e_local = 1:3
                v1, v2 = face_verts[e_local,f], face_verts[mod(e_local,3)+1,f]
                # find which global edge this corresponds to
                global_edge = findfirst(
                    x -> (x[1]==v1 && x[2]==v2) || (x[1]==v2 && x[2]==v1),
                    eachcol(edge_verts)
                )
                if global_edge !== nothing
                    faces[4:3+PD,f] .= edge_offsets[global_edge]:(edge_offsets[global_edge]+PD-2)
                end
            end
        end
    end

    # --- Step 3: face interior DOFs ---
    # Only for PD ≥ 3
    if PD ≥ 3
        nface_interior = (PD-1)*(PD-2) ÷ 2
        for f = 1:nfaces
            faces[4+PD:end,f] .= offset:(offset + nface_interior -1)
            offset += nface_interior
        end
    end

    return faces
end

function dof_coordinates(e::Tet{Lagrange, PD}) where PD
    coords = vertex_coordinates(e)

    offset = 4

    if PD < 2
        # return coords
        # do nothing to coords
    elseif PD == 2
        X = zeros(3, 10)

        # Vertices
        X[:, 1] .= (0, 0, 0)   # node 1
        X[:, 2] .= (1, 0, 0)   # node 2
        X[:, 3] .= (0, 1, 0)   # node 3
        X[:, 4] .= (0, 0, 1)   # node 4
    
        # Mid-edge nodes (ExodusII ordering)
        X[:, 5] .= (1//2, 0,     0    )  # edge 1-2
        X[:, 6] .= (1//2, 1//2,  0    )  # edge 2-3
        X[:, 7] .= (0,    1//2,  0    )  # edge 3-1
        X[:, 8] .= (0,    0,     1//2 )  # edge 1-4
        X[:, 9] .= (1//2, 0,     1//2 )  # edge 2-4
        X[:, 10] .= (0,    1//2,  1//2 )  # edge 3-4
    
        # return X
        coords = X
    else

        # edge DOFs
        if PD ≥ 2
            for (v1, v2) in eachcol(edge_vertices(e))
                for i in 1:PD - 1
                    t = i / PD
                    new_coord = (1.0 - t) * coords[:, v1] + t * coords[:, v2]
                    coords = hcat(coords, new_coord)
                    offset += 1
                end
            end
        end

        # face DOFs
        if PD ≥ 3
            for (v1, v2, v3) in face_vertices(e)
                for i in 1:PD - 1
                    for j in 1:PD - 1 - i
                        t1 = i / PD
                        t2 = j / PD
                        t3 = 1.0 - t1 - t2
                        new_coord = t1 * coords[:,v1] + t2 * coords[:,v2] + t3 * coords[:,v3]
                        coords = hcat(coords, new_coord)
                        offset += 1
                    end
                end
            end
        end

        # interior DOFs
        if PD ≥ 4
            for i in 1:PD - 3
                for j in 1:PD - 2 - i
                    for k in 1:PD - 1 - i - j
                        t1 = i / PD
                        t2 = j / PD
                        t3 = k / PD
                        t4 = 1.0 - t1 - t2 - t3
                        new_coord = t1 * coords[:,1] + t2 * coords[:,2] + t3 * coords[:,3] + t4 * coords[:,4]
                        coords = hcat(coords, new_coord)
                        offset += 1
                    end
                end
            end
        end
    end

    return coords
end

function interior_dofs(::Tet{Lagrange, PD}) where PD
    if PD < 4
        return Int[]
    else
        n_vertex = 4
        n_edge = 6 * max(PD - 1, 0)
        n_face = 4 * (PD - 1) * (PD - 2) ÷ 2
        start = n_vertex + n_edge + n_face + 1
        stop  = num_cell_dofs(Tet{Lagrange, PD}())
        return collect(start:stop)
    end
end

num_cell_dofs(::Tet{Lagrange, PD}) where PD = (PD + 1) * (PD + 2) * (PD + 3) ÷ 6
num_interior_dofs(::Tet{Lagrange, PD}) where PD = PD < 4 ? 0 : (PD - 1) * (PD - 2) * (PD - 3) ÷ 6

function cell_quadrature_points_and_weights(e::AbstractTet, q_rule::GaussLobattoLegendre)
    if cell_quadrature_degree(q_rule) == 1
        ξs = Matrix{Float64}(undef, 3, 1)
        ξs[:, 1] = [1. / 4., 1. / 4., 1. / 4.]
        ws = [1. / 6.]
    elseif cell_quadrature_degree(q_rule) == 2
        ξs = Matrix{Float64}(undef, 3, 5)
        ξs[:, 1] = [1. / 4., 1. / 4., 1. / 4.]
        ξs[:, 2] = [1. / 6., 1. / 6., 1. / 6.]
        ξs[:, 3] = [1. / 6., 1. / 6., 1. / 2.]
        ξs[:, 4] = [1. / 6., 1. / 2., 1. / 6.]
        ξs[:, 5] = [1. / 2., 1. / 6., 1. / 6.]

        #
        ws = [
            -2. / 15.
            3. / 40.
            3. / 40.
            3. / 40.
            3. / 40.
        ]
    else
        @assert false "Quadrature 1 through 2 currently supported."
    end
    return ξs, ws
end
function surface_quadrature_points_and_weights(e::AbstractTet, q_rule::GaussLobattoLegendre)
    ξs, ws = cell_quadrature_points_and_weights(boundary_element(e, 0), q_rule)

    ξ_return = zeros(3, length(ws), 4)
    w_return = zeros(length(ws), 4)

    ξ_return[1, :, 1] .= ξs[1, :]
    ξ_return[2, :, 1] .= 0.
    ξ_return[3, :, 1] .= ξs[2, :]
    #
    ξ_return[1, :, 2] .= ξs[1, :]
    ξ_return[2, :, 2] .= ξs[2, :]
    ξ_return[3, :, 2] .= 1. .- ξs[1, :] .- ξs[2, :]
    #
    ξ_return[1, :, 3] .= 0.
    ξ_return[2, :, 3] .= ξs[1, :]
    ξ_return[3, :, 3] .= ξs[2, :]
    #
    ξ_return[1, :, 4] .= ξs[1, :]
    ξ_return[2, :, 4] .= ξs[2, :]
    ξ_return[3, :, 4] .= 0.

    for n in 1:4
        w_return[:, n] .= ws
    end

    return ξ_return, w_return
end 
 
function shape_function_value(::Tet{Lagrange, 0}, _, _)
    return ones(1)
end

function shape_function_gradient(::Tet{Lagrange, 0}, _, _)
    return zeros(1, 3)
end

function shape_function_hessian(::Tet{Lagrange, 0}, _, _)
    return zeros(1, 3, 3)
end

function shape_function_value(::Tet{Lagrange, 1}, _, ξ)
    return [
        1. - ξ[1] - ξ[2] - ξ[3],
        ξ[1],
        ξ[2],
        ξ[3]
    ]
end

function shape_function_gradient(::Tet{Lagrange, 1}, _, ξ)
    grads = zeros(4, 3)

    grads[1, 1] = -1.
    grads[1, 2] = -1.
    grads[1, 3] = -1.
    #
    grads[2, 1] = 1.
    #
    grads[3, 2] = 1.
    #
    grads[4, 3] = 1.
    return grads
end

function shape_function_hessian(::Tet{Lagrange, 1}, _, ξ)
    return zeros(4, 3, 3)
end

function shape_function_value(::Tet{Lagrange, 2}, _, ξ)
    t0 = 1 - ξ[1] - ξ[2] - ξ[3]
    t1 = ξ[1]
    t2 = ξ[2]
    t3 = ξ[3]
    Ns = [
        t0 * (2 * t0 - 1),
        t1 * (2 * t1 - 1),
        t2 * (2 * t2 - 1),
        t3 * (2 * t3 - 1),
        4 * t0 * t1,
        4 * t1 * t2,
        4 * t2 * t0,
        4 * t0 * t3,
        4 * t1 * t3,
        4 * t2 * t3
    ]
    return Ns
end

function shape_function_gradient(::Tet{Lagrange, 2}, _, ξ)
    t0 = 1 - ξ[1] - ξ[2] - ξ[3]
    t1 = ξ[1]
    t2 = ξ[2]
    t3 = ξ[3]
    grads = zeros(10, 3)
    #
    grads[1, 1] = 1. - 4. * t0
    grads[2, 1] = 4. * t1 - 1
    grads[3, 1] = 0.
    grads[4, 1] = 0.
    grads[5, 1] = 4. * (t0 - t1)
    grads[6, 1] = 4. * t2
    grads[7, 1] = -4. * t2
    grads[8, 1] = -4. * t3
    grads[9, 1] = 4. * t3
    grads[10, 1] = 0.
    #
    grads[1, 2] = 1. - 4. * t0
    grads[2, 2] = 0.
    grads[3, 2] = 4. * t2 - 1.
    grads[4, 2] = 0.
    grads[5, 2] = -4. * t1
    grads[6, 2] = 4. * t1
    grads[7, 2] = 4. * (t0 - t2)
    grads[8, 2] = -4. * t3
    grads[9, 2] = 0.
    grads[10, 2] = 4. * t3
    #
    grads[1, 3] = 1. - 4. * t0
    grads[2, 3] = 0.
    grads[3, 3] = 0.
    grads[4, 3] = 4. * t3 - 1.
    grads[5, 3] = -4. * t1
    grads[6, 3] = 0.
    grads[7, 3] = -4. * t2
    grads[8, 3] = 4. * (t0 - t3)
    grads[9, 3] = 4. * t1
    grads[10, 3] = 4. * t2

    return grads
end

function shape_function_hessian(::Tet{Lagrange, 2}, _, _)
    hess = zeros(10, 3, 3)
    hess[:, 1, 1] .= [4.,  4.,  0.,  0., -8.,  0.,  0.,  0.,  0.,  0.]
    hess[:, 1, 2] .= [4.,  0.,  0.,  0., -4.,  4., -4.,  0.,  0.,  0.]
    hess[:, 1, 3] .= [4.,  0.,  0.,  0., -4.,  0.,  0., -4.,  4.,  0.]
    hess[:, 2, 1] .= [4.,  0.,  0.,  0., -4.,  4., -4.,  0.,  0.,  0.]
    hess[:, 2, 2] .= [4.,  0.,  4.,  0.,  0.,  0., -8.,  0.,  0.,  0.]
    hess[:, 2, 3] .= [4.,  0.,  0.,  0.,  0.,  0., -4., -4.,  0.,  4.]
    hess[:, 3, 1] .= [4.,  0.,  0.,  0., -4.,  0.,  0., -4.,  4.,  0.]
    hess[:, 3, 2] .= [4.,  0.,  0.,  0.,  0.,  0., -4., -4.,  0.,  4.]
    hess[:, 3, 3] .= [4.,  0.,  0.,  4.,  0.,  0.,  0., -8.,  0.,  0.]
    return hess
end

function shape_function_value(e::Tet{Lagrange, PD}, _, ξ) where PD
    # barycentric coordinates
    λ1 = 1 - ξ[1] - ξ[2] - ξ[3]
    λ2 = ξ[1]
    λ3 = ξ[2]
    λ4 = ξ[3]

    N = Vector{eltype(ξ)}(undef, num_cell_dofs(e))
    offset = 0

    # -------------------------
    # vertices
    # -------------------------
    offset += 1; N[offset] = λ1^PD
    offset += 1; N[offset] = λ2^PD
    offset += 1; N[offset] = λ3^PD
    offset += 1; N[offset] = λ4^PD

    # -------------------------
    # edges
    # -------------------------
    # edges in your canonical order: (V1,V2),(V2,V3),(V3,V1),(V1,V4),(V2,V4),(V3,V4)
    edges = [
        (λ1, λ2),
        (λ2, λ3),
        (λ3, λ1),
        (λ1, λ4),
        (λ2, λ4),
        (λ3, λ4)
    ]

    for (λa, λb) in edges
        for i in 1:PD-1
            offset += 1
            N[offset] = binomial(PD, i) * λa^(PD-i) * λb^i
        end
    end

    # -------------------------
    # faces
    # -------------------------
    # faces in your canonical order (face_vertices convention)
    faces = [
        (λ1, λ3, λ2),  # face 1
        (λ1, λ2, λ4),  # face 2
        (λ2, λ3, λ4),  # face 3
        (λ1, λ4, λ3)   # face 4
    ]

    for (λa, λb, λc) in faces
        for i in 1:PD-1
            for j in 1:PD-1-i
                k = PD - i - j
                offset += 1
                N[offset] = binomial(PD, i) *
                            binomial(PD - i, j) *
                            λa^i * λb^j * λc^k
            end
        end
    end

    # -------------------------
    # interior
    # -------------------------
    for i in 1:PD-3
        for j in 1:PD-2-i
            for k in 1:PD-1-i-j
                l = PD - i - j - k
                offset += 1
                N[offset] = binomial(PD, i) *
                            binomial(PD - i, j) *
                            binomial(PD - i - j, k) *
                            λ1^i * λ2^j * λ3^k * λ4^l
            end
        end
    end

    return N
end

function shape_function_gradient(e::Tet{Lagrange, PD}, _, ξ) where PD
    λ1 = 1 - ξ[1] - ξ[2] - ξ[3]
    λ2 = ξ[1]
    λ3 = ξ[2]
    λ4 = ξ[3]

    # barycentric gradients
    gλ1 = SVector(-1.0, -1.0, -1.0)
    gλ2 = SVector( 1.0,  0.0,  0.0)
    gλ3 = SVector( 0.0,  1.0,  0.0)
    gλ4 = SVector( 0.0,  0.0,  1.0)

    ndofs = num_cell_dofs(e)
    # dN = MVector{ndofs, SVector{3, eltype(ξ)}}(undef)  # 3D gradient
    dN = zeros(eltype(ξ), ndofs, 3)

    offset = 0

    # -------------------------
    # vertices
    # -------------------------
    λ = (λ1, λ2, λ3, λ4)
    gλ = (gλ1, gλ2, gλ3, gλ4)

    for i in 1:4
        offset += 1
        dN[offset, :] = PD * λ[i]^(PD - 1) * gλ[i]
    end

    # -------------------------
    # edges
    # -------------------------
    edges = edge_vertices(e)

    for (ia, ib) in eachcol(edges)
        λa, λb = λ[ia], λ[ib]
        ga, gb = gλ[ia], gλ[ib]
        for i in 1:PD-1
            offset += 1
            dN[offset] = binomial(PD, i) * (
                (PD-i) * λa^(PD-i-1) * λb^i * ga +
                i      * λa^(PD-i)   * λb^(i-1) * gb
            )
        end
    end

    # -------------------------
    # faces
    # -------------------------
    faces = [
        (1,3,2),
        (1,2,4),
        (2,3,4),
        (1,4,3)
    ]

    for (ia, ib, ic) in faces
        λa, λb, λc = λ[ia], λ[ib], λ[ic]
        ga, gb, gc = gλ[ia], gλ[ib], gλ[ic]
        for i in 1:PD-1
            for j in 1:PD-1-i
                k = PD - i - j
                offset += 1
                coeff = binomial(PD, i) * binomial(PD-i, j)
                dN[offset] = coeff * (
                    i * λa^(i-1) * λb^j * λc^k * ga +
                    j * λa^i     * λb^(j-1) * λc^k * gb +
                    k * λa^i     * λb^j     * λc^(k-1) * gc
                )
            end
        end
    end

    # -------------------------
    # interior DOFs
    # -------------------------
    for i in 1:PD-3
        for j in 1:PD-2-i
            for k in 1:PD-1-i-j
                l = PD - i - j - k
                offset += 1
                λa, λb, λc, λd = λ1, λ2, λ3, λ4
                ga, gb, gc, gd = gλ1, gλ2, gλ3, gλ4
                coeff = binomial(PD, i) * binomial(PD-i, j) * binomial(PD-i-j, k)
                dN[offset] = coeff * (
                    i * λa^(i-1) * λb^j * λc^k * λd^l * ga +
                    j * λa^i     * λb^(j-1) * λc^k * λd^l * gb +
                    k * λa^i     * λb^j     * λc^(k-1) * λd^l * gc +
                    l * λa^i     * λb^j     * λc^k     * λd^(l-1) * gd
                )
            end
        end
    end

    return dN
end
