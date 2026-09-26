# The surface Jacobian of MappedH1OrL2SurfaceInterpolants must turn the
# reference face measure into the physical face area: the quadrature sum of
# JxW over a face equals the area of that face.  For triangular faces the
# Jacobian was twice the area (|t1 × t2| is twice the triangle area and was
# divided by the reference area 1/2 once more), so every traction on a
# tetrahedron face was doubled.

function _face_area(X, ids)
    # polygon area from the vertices in order (triangle or planar quadrilateral)
    n = length(ids)
    a = zero(eltype(X)) * zeros(3)
    for k in 1:n
        i, j = ids[k], ids[mod1(k + 1, n)]
        a .+= cross(X[:, i], X[:, j])
    end
    return norm(a) / 2
end

# node coordinates of the reference element, one column per dof
function _node_coordinates(el::Tet{Lagrange, 1})
    return Float64.(vertex_coordinates(el))
end
function _node_coordinates(el::Tet{Lagrange, 2})
    V = Float64.(vertex_coordinates(el))
    M = reduce(hcat, [(V[:, a] .+ V[:, b]) ./ 2 for (a, b) in eachcol(edge_vertices(el))])
    return hcat(V, M)
end
function _node_coordinates(el::Tet{EnrichedLagrange, 2})
    return Float64.(reduce(hcat, ReferenceFiniteElements._nodes(el)))
end
function _node_coordinates(el::Hex{Lagrange, 1})
    return Float64.(vertex_coordinates(el))
end

function test_surface_jacobian(el, X)
    re = ReferenceFE(el, GaussLegendre(2))
    fv = face_vertices(el)
    NN = num_cell_dofs(el)
    X_el = SMatrix{3, NN, Float64}(X)
    for f in 1:size(fv, 2)
        area = 0.0
        for q in 1:num_surface_quadrature_points(re)
            area += MappedH1OrL2SurfaceInterpolants(re, X_el, q, f).JxW
        end
        @test area ≈ _face_area(X, fv[:, f]) rtol = 1e-12
    end
end

@testset "Surface Jacobian equals the face area" begin
    # tetrahedra: the reference element and an affine image of it
    A = [1.3 0.2 -0.4; 0.1 0.9 0.3; -0.2 0.5 1.7]
    for el in (Tet{Lagrange, 1}(), Tet{Lagrange, 2}(), Tet{EnrichedLagrange, 2}())
        Xref = _node_coordinates(el)
        test_surface_jacobian(el, Xref)
        test_surface_jacobian(el, A * Xref)
    end
    # hexahedron: the reference cube and an affine image (parallelogram faces)
    el = Hex{Lagrange, 1}()
    Xref = _node_coordinates(el)
    test_surface_jacobian(el, Xref)
    test_surface_jacobian(el, A * Xref)
end

# The surface quadrature points of face f must lie on face f as
# boundary_dofs defines it: the shape functions of the face's own nodes then
# sum to one there.  On the hexahedron the points of faces 1, 3, 5 and 6 lay
# on the faces 5, 6, 1 and 3, so a traction on those faces was integrated
# with a fraction of its value.  The tabulated boundary normal must also be
# the geometric normal of that face.
@testset "Surface quadrature points lie on their faces" begin
    for el in (Tet{Lagrange, 1}(), Tet{Lagrange, 2}(), Tet{EnrichedLagrange, 2}(),
               Hex{Lagrange, 1}(), Quad{Lagrange, 1}(), Quad{Lagrange, 2}(),
               Tri{Lagrange, 1}(), Tri{Lagrange, 2}())
        re = ReferenceFE(el, GaussLegendre(2))
        for f in 1:num_boundaries(el)
            bd = boundary_dofs(re, f)
            for q in 1:num_surface_quadrature_points(re)
                N = ReferenceFiniteElements.surface_shape_function_value(re, q, f)
                @test sum(N[bd]) ≈ 1 atol = 1e-12
            end
        end
    end
    # normals of the three-dimensional elements from the face geometry
    for el in (Tet{Lagrange, 1}(), Hex{Lagrange, 1}())
        V = Float64.(vertex_coordinates(el))
        fv = face_vertices(el)
        for f in 1:size(fv, 2)
            a, b, c = V[:, fv[1, f]], V[:, fv[2, f]], V[:, fv[3, f]]
            n = cross(b - a, c - a); n /= norm(n)
            @test n ≈ boundary_normals(el)[:, f] atol = 1e-12
        end
    end
end
