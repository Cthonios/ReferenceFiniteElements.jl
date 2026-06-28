"""
$(TYPEDEF)
"""
abstract type AbstractQuadratureType{CD, SD} end
"""
$(TYPEDSIGNATURES)
"""
function cell_quadrature_points_and_weights end
"""
$(TYPEDSIGNATURES)
"""
function surface_quadrature_points_and_weights end
"""
$(TYPEDSIGNATURES)
"""
cell_quadrature_degree(::AbstractQuadratureType{CD, SD}) where {CD, SD} = CD
"""
$(TYPEDSIGNATURES)
"""
surface_quadrature_degree(::AbstractQuadratureType{CD, SD}) where {CD, SD} = SD

"""
$(TYPEDEF)
"""
struct GaussLegendre{CD, SD} <: AbstractQuadratureType{CD, SD}
  function GaussLegendre(degree::Int)
    return GaussLegendre{degree, degree}()
  end

  function GaussLegendre(cell_degree::Int, surf_degree::Int)
    return GaussLegendre{cell_degree, surf_degree}()
  end

  function GaussLegendre{CD, SD}() where {CD, SD}
    @assert isa(CD, Integer)
    @assert isa(SD, Integer)
    @assert CD > 0 "Cell quadrature degree must be greater than zero"
    @assert SD > 0 "Surface quadrature degree must be greater than zero"
    new{CD, SD}()
  end
end

########################################################################
# Edge implementation
########################################################################
function cell_quadrature_points_and_weights(e::AbstractEdge, q_rule::GaussLegendre)
    ξs, ws = gausslegendre(cell_quadrature_degree(q_rule))

    # if e.shifted
    if _is_shifted(e)
        ξs .= (ξs .+ 1.) ./ 2.
        ws .= ws ./ 2.
    end

    return reshape(ξs, 1, length(ξs)), ws
end

num_cell_quadrature_points(::AbstractEdge, ::Type{GaussLegendre{CD, SD}}) where {CD, SD} = CD

function surface_quadrature_points_and_weights(e::AbstractEdge, ::GaussLegendre)
    if _is_shifted(e)
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

########################################################################
# Hex implementation
########################################################################
function cell_quadrature_points_and_weights(e::AbstractHex, q_rule::GaussLegendre)
    ξs, ws = cell_quadrature_points_and_weights(boundary_element(boundary_element(e, 0), 0), q_rule)
    n = length(ws)
    ξ_return = Matrix{eltype(ξs)}(undef, 3, n * n * n)
    w_return = Vector{eltype(ξs)}(undef, n * n * n)
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

num_cell_quadrature_points(::AbstractHex, ::Type{GaussLegendre{CD, SD}}) where {CD, SD} = CD * CD * CD

function surface_quadrature_points_and_weights(e::AbstractHex, q_rule::GaussLegendre)
    ξs, ws = cell_quadrature_points_and_weights(boundary_element(e, 0), q_rule)

    ξ_return = zeros(3, length(ws), 6)
    w_return = zeros(length(ws), 6)

    ξ_return[1:2, :, 1] .= ξs
    ξ_return[3, :, 1]   .= -1.
    ξ_return[1, :, 2]   .= 1.
    ξ_return[2:3, :, 2] .= ξs
    ξ_return[1:2, :, 3] .= ξs
    ξ_return[3, :, 3]   .= 1.
    ξ_return[1, :, 4]   .= -1.
    ξ_return[2:3, :, 4] .= ξs
    ξ_return[1, :, 5]   .= ξs[1, :]
    ξ_return[2, :, 5]   .= -1.
    ξ_return[3, :, 5]   .= ξs[2, :]
    ξ_return[1, :, 6]   .= ξs[1, :]
    ξ_return[2, :, 6]   .= 1.
    ξ_return[3, :, 6]   .= ξs[2, :]

    for n in 1:6
        w_return[:, n] .= ws
    end
    return ξ_return, w_return
end

########################################################################
# Quad implementation
########################################################################
function cell_quadrature_points_and_weights(e::AbstractQuad, q_rule::GaussLegendre)
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

num_cell_quadrature_points(::AbstractQuad, ::Type{GaussLegendre{CD, SD}}) where {CD, SD} = CD * CD

function surface_quadrature_points_and_weights(e::AbstractQuad, q_rule::GaussLegendre)
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

########################################################################
# Tet implementation
########################################################################
function cell_quadrature_points_and_weights(::AbstractTet, q_rule::GaussLegendre)
    deg = cell_quadrature_degree(q_rule)
    if deg == 1
        # 1-point centroid rule (degree 1)
        ξs = Matrix{Float64}(undef, 3, 1)
        ξs[:, 1] = [1. / 4., 1. / 4., 1. / 4.]
        ws = [1. / 6.]
    elseif deg == 2
        # 4-point symmetric rule (degree 2)
        s = sqrt(5.0)
        a = (5. + 3. * s) / 20.
        b = (5. - s) / 20.
        ξs = Matrix{Float64}(undef, 3, 4)
        ξs[:, 1] = [b, b, b]
        ξs[:, 2] = [a, b, b]
        ξs[:, 3] = [b, a, b]
        ξs[:, 4] = [b, b, a]
        ws = [1. / 24., 1. / 24., 1. / 24., 1. / 24.]
    elseif deg == 3
        # 5-point rule (degree 3)
        ξs = Matrix{Float64}(undef, 3, 5)
        ξs[:, 1] = [1. / 4., 1. / 4., 1. / 4.]
        ξs[:, 2] = [1. / 6., 1. / 6., 1. / 6.]
        ξs[:, 3] = [1. / 6., 1. / 6., 1. / 2.]
        ξs[:, 4] = [1. / 6., 1. / 2., 1. / 6.]
        ξs[:, 5] = [1. / 2., 1. / 6., 1. / 6.]
        ws = [-2. / 15., 3. / 40., 3. / 40., 3. / 40., 3. / 40.]
    else
        @assert false "GaussLegendre degree 1 through 3 supported for Tet."
    end
    return ξs, ws
end

num_cell_quadrature_points(::AbstractTet, ::Type{GaussLegendre{1, SD}}) where SD = 1
num_cell_quadrature_points(::AbstractTet, ::Type{GaussLegendre{2, SD}}) where SD = 4
num_cell_quadrature_points(::AbstractTet, ::Type{GaussLegendre{3, SD}}) where SD = 5

function surface_quadrature_points_and_weights(e::AbstractTet, q_rule::GaussLegendre)
    return surface_quadrature_points_and_weights(e, GaussLobattoLegendre(cell_quadrature_degree(q_rule), surface_quadrature_degree(q_rule)))
end

########################################################################
# Tri implementation
########################################################################
function cell_quadrature_points_and_weights(::AbstractTri, q_rule::GaussLegendre)
    deg = cell_quadrature_degree(q_rule)
    if deg == 1
        # 1-point centroid rule (degree 1)
        ξs = Matrix{Float64}(undef, 2, 1)
        ξs[:, 1] = [1. / 3., 1. / 3.]
        ws = [0.5]
    elseif deg == 2
        # 3-point rule (degree 2)
        ξs = Matrix{Float64}(undef, 2, 3)
        ξs[:, 1] = [1. / 6., 1. / 6.]
        ξs[:, 2] = [4. / 6., 1. / 6.]
        ξs[:, 3] = [1. / 6., 4. / 6.]
        ws = [1. / 6., 1. / 6., 1. / 6.]
    elseif deg == 3
        # 4-point rule (degree 3): centroid + 3 edge midpoints
        ξs = Matrix{Float64}(undef, 2, 4)
        ξs[:, 1] = [1. / 3., 1. / 3.]
        ξs[:, 2] = [1. / 5., 3. / 5.]
        ξs[:, 3] = [3. / 5., 1. / 5.]
        ξs[:, 4] = [1. / 5., 1. / 5.]
        ws = [-27. / 96., 25. / 96., 25. / 96., 25. / 96.]
    else
        @assert false "GaussLegendre degree 1 through 3 supported for Tri."
    end
    return ξs, ws
end

num_cell_quadrature_points(::AbstractTri, ::Type{GaussLegendre{1, SD}}) where SD = 1
num_cell_quadrature_points(::AbstractTri, ::Type{GaussLegendre{2, SD}}) where SD = 3
num_cell_quadrature_points(::AbstractTri, ::Type{GaussLegendre{3, SD}}) where SD = 4

function surface_quadrature_points_and_weights(e::AbstractTri, q_rule::GaussLegendre)
    return surface_quadrature_points_and_weights(e, GaussLobattoLegendre(cell_quadrature_degree(q_rule), surface_quadrature_degree(q_rule)))
end

"""
$(TYPEDEF)
"""
struct GaussLobattoLegendre{CD, SD} <: AbstractQuadratureType{CD, SD}
  function GaussLobattoLegendre(degree::Int)
    return GaussLobattoLegendre{degree, degree}()
  end
  
  function GaussLobattoLegendre(cell_degree::Int, surf_degree::Int)
    return GaussLobattoLegendre{cell_degree, surf_degree}()
  end

  function GaussLobattoLegendre{CD, SD}() where {CD, SD}
    @assert isa(CD, Integer)
    @assert isa(SD, Integer)
    @assert CD > 0 "Cell quadrature degree must be greater than zero"
    @assert SD > 0 "Surface quadrature degree must be greater than zero"
    new{CD, SD}()
  end
end

########################################################################
# Edge implementation
########################################################################
function cell_quadrature_points_and_weights(e::AbstractEdge, q_rule::GaussLobattoLegendre)
    ξs, ws = gausslegendre(cell_quadrature_degree(q_rule))

    if _is_shifted(e)
        ξs .= (ξs .+ 1.) ./ 2.
        ws .= ws ./ 2.
    end

    return reshape(ξs, 1, length(ξs)), ws
end

num_cell_quadrature_points(::AbstractEdge, ::Type{GaussLobattoLegendre{CD, SD}}) where {CD, SD} = CD

function surface_quadrature_points_and_weights(e::AbstractEdge, ::GaussLobattoLegendre)
    if _is_shifted(e)
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

########################################################################
# Hex implementation
########################################################################
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

num_cell_quadrature_points(::AbstractHex, ::Type{GaussLobattoLegendre{CD, SD}}) where {CD, SD} = CD * CD * CD

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

########################################################################
# Quad implementation
########################################################################
function cell_quadrature_points_and_weights(e::AbstractQuad, q_rule::GaussLobattoLegendre)
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

num_cell_quadrature_points(::AbstractQuad, ::Type{GaussLobattoLegendre{CD, SD}}) where {CD, SD} = CD * CD

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

########################################################################
# Tet implementation
########################################################################
function cell_quadrature_points_and_weights(::AbstractTet, q_rule::GaussLobattoLegendre)
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

num_cell_quadrature_points(::AbstractTet, ::Type{GaussLobattoLegendre{1, SD}}) where SD = 1
num_cell_quadrature_points(::AbstractTet, ::Type{GaussLobattoLegendre{2, SD}}) where SD = 5

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

########################################################################
# Tri implementation
########################################################################
function cell_quadrature_points_and_weights(::AbstractTri, q_rule::GaussLobattoLegendre)
    if cell_quadrature_degree(q_rule) == 1
      ξs = Matrix{Float64}(undef, 2, 1)
      ξs[:, 1] = [1. / 3., 1. / 3.]
      ws = [0.5]
    elseif cell_quadrature_degree(q_rule) == 2
      ξs = Matrix{Float64}(undef, 2, 3)
      ξs[:, 1] = [2. / 3., 1. / 6.]
      ξs[:, 2] = [1. / 6., 2. / 3.]
      ξs[:, 3] = [1. / 6., 1. / 6.]
      ws = [1. / 6., 1. / 6., 1. / 6.]
    elseif cell_quadrature_degree(q_rule) <= 4
      ξs = Matrix{Float64}(undef, 2, 6)
      ξs[:, 1] = [1.081030181680700E-01, 4.459484909159650E-01]
      ξs[:, 2] = [4.459484909159650E-01, 1.081030181680700E-01]
      ξs[:, 3] = [4.459484909159650E-01, 4.459484909159650E-01]
      ξs[:, 4] = [8.168475729804590E-01, 9.157621350977100E-02]
      ξs[:, 5] = [9.157621350977100E-02, 8.168475729804590E-01]
      ξs[:, 6] = [9.157621350977100E-02, 9.157621350977100E-02]
  
      ws = [
        1.116907948390055E-01,
        1.116907948390055E-01,
        1.116907948390055E-01,
        5.497587182766100E-02,
        5.497587182766100E-02,
        5.497587182766100E-02
      ]
    elseif cell_quadrature_degree(q_rule) <= 5
      ξs = Matrix{Float64}(undef, 2, 7)
      ξs[:, 1] = [3.33333333333333E-01, 3.33333333333333E-01]
      ξs[:, 2] = [5.97158717897700E-02, 4.70142064105115E-01]
      ξs[:, 3] = [4.70142064105115E-01, 5.97158717897700E-02]
      ξs[:, 4] = [4.70142064105115E-01, 4.70142064105115E-01]
      ξs[:, 5] = [7.97426985353087E-01, 1.01286507323456E-01]
      ξs[:, 6] = [1.01286507323456E-01, 7.97426985353087E-01]
      ξs[:, 7] = [1.01286507323456E-01, 1.01286507323456E-01]
  
      ws = [
        1.12500000000000E-01,
        6.61970763942530E-02,
        6.61970763942530E-02,
        6.61970763942530E-02,
        6.29695902724135E-02,
        6.29695902724135E-02,
        6.29695902724135E-02
      ]
    elseif cell_quadrature_degree(q_rule) <= 6
      ξs = Matrix{Float64}(undef, 2, 12)
      ξs[:, 1]  = [5.01426509658179E-01, 2.49286745170910E-01]
      ξs[:, 2]  = [2.49286745170910E-01, 5.01426509658179E-01]
      ξs[:, 3]  = [2.49286745170910E-01, 2.49286745170910E-01]
      ξs[:, 4]  = [8.73821971016996E-01, 6.30890144915020E-02]
      ξs[:, 5]  = [6.30890144915020E-02, 8.73821971016996E-01]
      ξs[:, 6]  = [6.30890144915020E-02, 6.30890144915020E-02]
      ξs[:, 7]  = [5.31450498448170E-02, 3.10352451033784E-01]
      ξs[:, 8]  = [6.36502499121399E-01, 5.31450498448170E-02]
      ξs[:, 9]  = [3.10352451033784E-01, 6.36502499121399E-01]
      ξs[:, 10] = [5.31450498448170E-02, 6.36502499121399E-01]
      ξs[:, 11] = [6.36502499121399E-01, 3.10352451033784E-01]
      ξs[:, 12] = [3.10352451033784E-01, 5.31450498448170E-02]
  
      ws = [
        5.83931378631895E-02,
        5.83931378631895E-02,
        5.83931378631895E-02,
        2.54224531851035E-02,
        2.54224531851035E-02,
        2.54224531851035E-02,
        4.14255378091870E-02,
        4.14255378091870E-02,
        4.14255378091870E-02,
        4.14255378091870E-02,
        4.14255378091870E-02,
        4.14255378091870E-02
      ]
    else
        @assert false "Quadrature degree 1 through 6 currently supported."
    end

    return ξs, ws
end

num_cell_quadrature_points(::AbstractTri, ::Type{GaussLobattoLegendre{1, SD}}) where SD = 1
num_cell_quadrature_points(::AbstractTri, ::Type{GaussLobattoLegendre{2, SD}}) where SD = 3
num_cell_quadrature_points(::AbstractTri, ::Type{GaussLobattoLegendre{3, SD}}) where SD = 6
num_cell_quadrature_points(::AbstractTri, ::Type{GaussLobattoLegendre{4, SD}}) where SD = 6
num_cell_quadrature_points(::AbstractTri, ::Type{GaussLobattoLegendre{5, SD}}) where SD = 7
num_cell_quadrature_points(::AbstractTri, ::Type{GaussLobattoLegendre{6, SD}}) where SD = 12

function surface_quadrature_points_and_weights(e::AbstractTri, q_rule::GaussLobattoLegendre)
    ξs, ws = cell_quadrature_points_and_weights(boundary_element(e, 0), q_rule)

    ξ_return = zeros(2, length(ws), 3)
    w_return = zeros(length(ws), 3)

    ξ_return[1, :, 1] .= ξs[1, :]
    ξ_return[2, :, 1] .= -1.
    ξ_return[1, :, 2] .= ξs[1, :]
    ξ_return[2, :, 2] .= 1. .- ξs[1, :]
    ξ_return[1, :, 3] .= -1.
    ξ_return[2, :, 3] .= ξs[1, :]

    for n in 1:3
        w_return[:, n] .= ws
    end
    return ξ_return, w_return
end
