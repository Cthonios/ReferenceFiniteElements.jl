# Lagrange spaces enriched with bubbles in a nodal basis: the seven-node
# triangle and the fifteen-node tetrahedron.

function _fd_gradient(el, ξ; h = 1e-6)
    d = length(ξ)
    g = zeros(num_cell_dofs(el), d)
    for j in 1:d
        e = zeros(d); e[j] = h
        g[:, j] = (ReferenceFiniteElements.shape_function_value(el, ξ .+ e) .-
                   ReferenceFiniteElements.shape_function_value(el, ξ .- e)) ./ (2h)
    end
    return g
end

function _fd_hessian(el, ξ; h = 1e-5)
    d = length(ξ)
    H = zeros(num_cell_dofs(el), d, d)
    for j in 1:d
        e = zeros(d); e[j] = h
        H[:, :, j] = (ReferenceFiniteElements.shape_function_gradient(el, ξ .+ e) .-
                      ReferenceFiniteElements.shape_function_gradient(el, ξ .- e)) ./ (2h)
    end
    return H
end

"Random points in the reference simplex of dimension d."
function _random_simplex_points(d, n)
    pts = Vector{Vector{Float64}}()
    while length(pts) < n
        ξ = rand(d)
        sum(ξ) < 1 && push!(pts, ξ)
    end
    return pts
end

"Random points on face f of the reference tetrahedron (RFE face_vertices order)."
function _random_face_points(f, n)
    verts = ([0., 0, 0], [1., 0, 0], [0., 1, 0], [0., 0, 1])
    a, b, c = face_vertices(Tet{Lagrange, 2}())[:, f]
    pts = Vector{Vector{Float64}}()
    for _ in 1:n
        s, t = rand(), rand()
        s + t > 1 && ((s, t) = (1 - s, 1 - t))
        push!(pts, (1 - s - t) .* verts[a] .+ s .* verts[b] .+ t .* verts[c])
    end
    return pts
end

function test_enriched_element(el, quadratic, d)
    n = num_cell_dofs(el)
    nodes = ReferenceFiniteElements._nodes(el)
    @test length(nodes) == n

    @testset "nodal basis: unit value at own node, zero at the others" begin
        for (j, ξ) in enumerate(nodes)
            N = ReferenceFiniteElements.shape_function_value(el, Float64.(ξ))
            for i in 1:n
                @test isapprox(N[i], i == j ? 1.0 : 0.0; atol = 1e-13)
            end
        end
    end

    @testset "partition of unity and its derivatives" begin
        for ξ in _random_simplex_points(d, 20)
            @test sum(ReferenceFiniteElements.shape_function_value(el, ξ)) ≈ 1.0
            @test all(abs.(sum(ReferenceFiniteElements.shape_function_gradient(el, ξ); dims = 1)) .< 1e-12)
            @test all(abs.(sum(ReferenceFiniteElements.shape_function_hessian(el, ξ); dims = 1)) .< 1e-11)
        end
    end

    @testset "gradient and hessian agree with finite differences" begin
        for ξ in _random_simplex_points(d, 10)
            g = ReferenceFiniteElements.shape_function_gradient(el, ξ)
            @test maximum(abs, g - _fd_gradient(el, ξ)) < 1e-8
            H = ReferenceFiniteElements.shape_function_hessian(el, ξ)
            @test maximum(abs, H - _fd_hessian(el, ξ)) < 1e-7
        end
    end

    @testset "the space contains the quadratics" begin
        # Interpolating a quadratic at the nodes reproduces it everywhere,
        # which is what makes the element at least as accurate as P2.
        coeffs = [quadratic(Float64.(ξ)) for ξ in nodes]
        for ξ in _random_simplex_points(d, 20)
            N = ReferenceFiniteElements.shape_function_value(el, ξ)
            @test dot(coeffs, N) ≈ quadratic(ξ) atol = 1e-12
        end
    end
end

function test_enriched_tet()
    el = Tet{EnrichedLagrange, 2}()
    quadratic(ξ) = 1 + 2ξ[1] - ξ[2] + 0.5ξ[3] + ξ[1]^2 - 3ξ[2]*ξ[3] + 0.25ξ[3]^2 + ξ[1]*ξ[2]
    test_enriched_element(el, quadratic, 3)

    @testset "dof interface" begin
        @test num_cell_dofs(el) == 15
        @test num_interior_dofs(el) == 1
        @test interior_dofs(el) == [15]
        bd = boundary_dofs(el)
        @test size(bd) == (7, 4)
        for f in 1:4
            @test num_dofs_on_boundary(el, f) == 7
            @test bd[7, f] == 10 + f
        end
        @test boundary_element(el, 1) == Tri{EnrichedLagrange, 2}()
        @test boundary_element(boundary_element(el, 1), 1) == Edge{Lagrange, 2}(; shifted = true)
    end

    @testset "trace: the functions of a face are the ones nonzero on it" begin
        bd = boundary_dofs(el)
        for f in 1:4
            on_face = Set(bd[:, f])
            for ξ in _random_face_points(f, 10)
                N = ReferenceFiniteElements.shape_function_value(el, ξ)
                for i in 1:15
                    i in on_face || @test abs(N[i]) < 1e-13
                end
            end
        end
    end

    @testset "the face nodes lie at the face centroids, in face_vertices order" begin
        verts = ([0., 0, 0], [1., 0, 0], [0., 1, 0], [0., 0, 1])
        fv = face_vertices(el)
        nodes = ReferenceFiniteElements._nodes(el)
        for f in 1:4
            @test Float64.(nodes[10 + f]) ≈ (verts[fv[1, f]] .+ verts[fv[2, f]] .+ verts[fv[3, f]]) ./ 3
        end
        @test Float64.(nodes[15]) ≈ [0.25, 0.25, 0.25]
    end

    @testset "degree-5 rule: exact through degree 5, not at 6" begin
        re = ReferenceFE(el, GaussLegendre(5))
        @test num_cell_quadrature_points(re) == 14
        worst = 0.0
        for deg in 0:5, a in 0:deg, b in 0:(deg - a)
            c = deg - a - b
            exact = factorial(a) * factorial(b) * factorial(c) / factorial(deg + 3)
            num = sum(cell_quadrature_weight(re, q) * prod(cell_quadrature_point(re, q) .^ (a, b, c))
                      for q in 1:14)
            worst = max(worst, abs(num - exact) / exact)
        end
        @test worst < 1e-14
        exact6 = factorial(6) / factorial(9)
        num6 = sum(cell_quadrature_weight(re, q) * cell_quadrature_point(re, q)[1]^6 for q in 1:14)
        @test abs(num6 - exact6) / exact6 > 1e-3
        for q in 1:14
            @test cell_quadrature_weight(re, q) > 0
        end
        @test ReferenceFE(el, GaussLegendre(4)) isa ReferenceFE
    end
end

function test_enriched_tri()
    el = Tri{EnrichedLagrange, 2}()
    quadratic(ξ) = 1 - ξ[1] + 2ξ[2] + ξ[1]^2 - ξ[1]*ξ[2] + 0.5ξ[2]^2
    test_enriched_element(el, quadratic, 2)
    @testset "dof interface" begin
        @test num_cell_dofs(el) == 7
        @test interior_dofs(el) == [7]
        @test boundary_dofs(el) == boundary_dofs(Tri{Lagrange, 2}())
        @test num_dofs_on_boundary(el, 1) == 3
        @test Float64.(ReferenceFiniteElements._nodes(el)[7]) ≈ [1/3, 1/3]
    end
    @testset "the centroid function vanishes on the edges" begin
        for s in rand(10)
            for ξ in ([s, 0.], [1 - s, s], [0., s])
                @test abs(ReferenceFiniteElements.shape_function_value(el, ξ)[7]) < 1e-13
            end
        end
    end
end

"For every element, the dofs listed for face f are exactly the functions that
are nonzero on the face_vertices face f; the surface quadrature points and the
boundary normals follow face_vertices, so a mismatch sends side-set data to
the wrong face."
function test_face_dofs_follow_face_vertices(el)
    verts = ([0., 0, 0], [1., 0, 0], [0., 1, 0], [0., 0, 1])
    fv = face_vertices(el)
    bd = boundary_dofs(el)
    n = num_cell_dofs(el)
    for f in 1:size(fv, 2)
        a, b, c = fv[:, f]
        on_face = Set(bd[:, f])
        nonzero = falses(n)
        for _ in 1:20
            s, t = rand(), rand()
            s + t > 1 && ((s, t) = (1 - s, 1 - t))
            ξ = (1 - s - t) .* verts[a] .+ s .* verts[b] .+ t .* verts[c]
            N = ReferenceFiniteElements.shape_function_value(el, ξ)
            nonzero .|= abs.(N) .> 1e-12
        end
        @test Set(findall(nonzero)) == on_face
    end
end

function test_tri_hessian_against_finite_differences()
    for p in 1:4
        el = Tri{Lagrange, p}()
        for ξ in _random_simplex_points(2, 5)
            H = ReferenceFiniteElements.shape_function_hessian(el, ξ)
            @test maximum(abs, H - _fd_hessian(el, ξ)) < 1e-7
        end
    end
end
