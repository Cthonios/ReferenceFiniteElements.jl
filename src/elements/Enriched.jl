########################################################################
# Lagrange spaces enriched with bubbles, in a nodal basis
########################################################################
#
# The space is the Lagrange space of degree 2 together with bubble functions,
# products of barycentric coordinates that vanish on a prescribed part of the
# element boundary.  Each element carries the space in a nodal basis: the
# function attached to a node has unit value there and vanishes at the other
# nodes, so that a degree-of-freedom manager keyed by node handles the element
# like any Lagrange element.  The nodal basis is obtained from the
# hierarchical one (Lagrange functions, then bubbles) through the constant
# matrix A = V^{-T}, with V the values of the hierarchical functions at the
# nodes; V has rational entries, so A is computed exactly and stored in
# floating point.
#
# Node ordering.
#   Tri,  degree 2, 7 nodes:  1-3 vertices, 4-6 edge midpoints in the order of
#         edge_vertices (1-2, 2-3, 3-1), 7 centroid.
#   Tet,  degree 2, 15 nodes: 1-4 vertices, 5-10 edge midpoints in the order
#         of edge_vertices (1-2, 2-3, 3-1, 1-4, 2-4, 3-4), 11-14 face
#         centroids in the order of face_vertices (1-2-4, 2-3-4, 1-4-3,
#         1-3-2), 15 centroid.
# The bubble of face f is 27 λa λb λc over the vertices of that face, with
# unit value at its centroid; the interior bubble of the tetrahedron is
# 256 λ1 λ2 λ3 λ4 and that of the triangle 27 λ1 λ2 λ3, each with unit
# value at the centroid.

# Value, gradient and hessian of c * prod(λ[m] for m in idx), with λ the
# barycentric coordinates and gλ their (constant) gradients.
function _bubble(idx, c, λ, gλ, ::Val{ND}) where ND
    T = promote_type(eltype(λ), Float64)
    v = c * prod(λ[m] for m in idx)
    g = zeros(T, ND)
    for m in idx
        f = c
        for n in idx
            n == m || (f *= λ[n])
        end
        g .+= f .* gλ[m]
    end
    h = zeros(T, ND, ND)
    for m in idx, n in idx
        m == n && continue
        f = c
        for k in idx
            (k == m || k == n) || (f *= λ[k])
        end
        for i in 1:ND, j in 1:ND
            h[i, j] += f * gλ[m][i] * gλ[n][j]
        end
    end
    return v, g, h
end

# ---------------------------------------------------------------------------
# Triangle: P2 plus the interior bubble, seven nodes
# ---------------------------------------------------------------------------
const _TRI7_BUBBLES = ((1, 2, 3),)

function _hierarchical(e::Tri{EnrichedLagrange, 2}, ξ)
    base = Tri{Lagrange, 2}()
    λ = (1 - ξ[1] - ξ[2], ξ[1], ξ[2])
    gλ = (SVector(-1.0, -1.0), SVector(1.0, 0.0), SVector(0.0, 1.0))
    N = collect(shape_function_value(base, ξ))
    dN = Matrix(shape_function_gradient(base, ξ))
    d2N = Array(shape_function_hessian(base, ξ))
    for idx in _TRI7_BUBBLES
        v, g, h = _bubble(idx, 27, λ, gλ, Val(2))
        push!(N, v)
        dN = vcat(dN, reshape(g, 1, 2))
        d2N = cat(d2N, reshape(h, 1, 2, 2); dims = 1)
    end
    return N, dN, d2N
end

_nodes(::Tri{EnrichedLagrange, 2}) = (
    [0, 0], [1, 0], [0, 1],
    [1//2, 0], [1//2, 1//2], [0, 1//2],
    [1//3, 1//3],
)

# ---------------------------------------------------------------------------
# Tetrahedron: P2 plus four face bubbles and the interior bubble, fifteen nodes
# ---------------------------------------------------------------------------
# face bubbles in the order of face_vertices, then the interior bubble
const _TET15_BUBBLES = ((1, 2, 4), (2, 3, 4), (1, 4, 3), (1, 3, 2), (1, 2, 3, 4))

function _hierarchical(e::Tet{EnrichedLagrange, 2}, ξ)
    base = Tet{Lagrange, 2}()
    λ = (1 - ξ[1] - ξ[2] - ξ[3], ξ[1], ξ[2], ξ[3])
    gλ = (SVector(-1.0, -1.0, -1.0), SVector(1.0, 0.0, 0.0),
          SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
    N = collect(shape_function_value(base, ξ))
    dN = Matrix(shape_function_gradient(base, ξ))
    d2N = Array(shape_function_hessian(base, ξ))
    for idx in _TET15_BUBBLES
        c = length(idx) == 3 ? 27 : 256
        v, g, h = _bubble(idx, c, λ, gλ, Val(3))
        push!(N, v)
        dN = vcat(dN, reshape(g, 1, 3))
        d2N = cat(d2N, reshape(h, 1, 3, 3); dims = 1)
    end
    return N, dN, d2N
end

function _nodes(::Tet{EnrichedLagrange, 2})
    v = ([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1])
    pts = Vector{Vector{Rational{Int}}}([Rational{Int}.(x) for x in v])
    for (a, b) in eachcol(edge_vertices(Tet{Lagrange, 2}()))
        push!(pts, (v[a] .+ v[b]) .// 2)
    end
    for (a, b, c) in eachcol(face_vertices(Tet{Lagrange, 2}()))
        push!(pts, (v[a] .+ v[b] .+ v[c]) .// 3)
    end
    push!(pts, [1//4, 1//4, 1//4])
    return pts
end

# ---------------------------------------------------------------------------
# The nodal transform, computed exactly once per element type
# ---------------------------------------------------------------------------
function _nodal_transform(e)
    pts = _nodes(e)
    n = length(pts)
    V = Matrix{Rational{BigInt}}(undef, n, n)
    for (j, ξ) in enumerate(pts)
        V[j, :] = _hierarchical(e, Rational{BigInt}.(ξ))[1]
    end
    A = inv(V)'
    return SMatrix{n, n, Float64}(Float64.(A))
end

const _TRI7_A  = _nodal_transform(Tri{EnrichedLagrange, 2}())
const _TET15_A = _nodal_transform(Tet{EnrichedLagrange, 2}())

_transform(::Tri{EnrichedLagrange, 2}) = _TRI7_A
_transform(::Tet{EnrichedLagrange, 2}) = _TET15_A

const _EnrichedElement = Union{Tri{EnrichedLagrange, 2}, Tet{EnrichedLagrange, 2}}

function shape_function_value(e::_EnrichedElement, ξ)
    N, _, _ = _hierarchical(e, ξ)
    return _transform(e) * N
end

function shape_function_gradient(e::_EnrichedElement, ξ)
    _, dN, _ = _hierarchical(e, ξ)
    return _transform(e) * dN
end

function shape_function_hessian(e::_EnrichedElement, ξ)
    _, _, d2N = _hierarchical(e, ξ)
    A = _transform(e)
    n, d = size(d2N, 1), size(d2N, 2)
    out = zeros(eltype(d2N), n, d, d)
    for i in 1:d, j in 1:d
        out[:, i, j] = A * d2N[:, i, j]
    end
    return out
end

# ---------------------------------------------------------------------------
# Degrees of freedom
# ---------------------------------------------------------------------------
num_cell_dofs(::Tri{EnrichedLagrange, 2}) = 7
num_interior_dofs(::Tri{EnrichedLagrange, 2}) = 1
interior_dofs(::Tri{EnrichedLagrange, 2}) = [7]
boundary_dofs(::Tri{EnrichedLagrange, 2}) = boundary_dofs(Tri{Lagrange, 2}())
num_dofs_on_boundary(::Tri{EnrichedLagrange, 2}, ::Int) = 3

num_cell_dofs(::Tet{EnrichedLagrange, 2}) = 15
num_interior_dofs(::Tet{EnrichedLagrange, 2}) = 1
interior_dofs(::Tet{EnrichedLagrange, 2}) = [15]
function boundary_dofs(::Tet{EnrichedLagrange, 2})
    p2 = boundary_dofs(Tet{Lagrange, 2}())          # 6 x 4
    return vcat(p2, reshape(collect(11:14), 1, 4))   # face node of each face
end
