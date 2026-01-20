"""
$(TYPEDEF)
"""
abstract type AbstractInterpolants end
abstract type AbstractDynamicInterpolants <: AbstractInterpolants end # TODO implment
abstract type AbstractStaticInterpolants <: AbstractInterpolants end
dimension(interps::AbstractStaticInterpolants) = size(interps.ξ, 1)
num_quadrature_points(interps::AbstractVector{T}) where T <: AbstractStaticInterpolants = length(interps)
num_quadrature_points(interps::AbstractMatrix{T}) where T <: AbstractStaticInterpolants = size(interps, 1)
quadrature_point(interps::AbstractVector{T}, q::Int) where T <: AbstractStaticInterpolants = interps[q].ξ
quadrature_point(interps::AbstractMatrix{T}, q::Int, f::Int) where T <: AbstractStaticInterpolants = interps[q, f].ξ
quadrature_weight(interps::AbstractVector{T}, q::Int) where T <: AbstractStaticInterpolants = interps[q].w
quadrature_weight(interps::AbstractMatrix{T}, q::Int, f::Int) where T <: AbstractStaticInterpolants = interps[q, f].w
shape_function_value(interps::AbstractVector{T}, q::Int) where T <: AbstractStaticInterpolants = interps[q].N
shape_function_value(interps::AbstractMatrix{T}, q::Int, f::Int) where T <: AbstractStaticInterpolants = interps[q, f].N

# abstract type AbstractH1OrL2Interpolants <: AbstractInterpolants end
# abstract type AbstractHdivInterpolants <: AbstractInterpolants end

# struct H1OrL2Interpolants{
#     RT <: Number,
#     R1 <: Union{<:AbstractArray{RT, 1}, <:AbstractArray{RT, 2}},
#     R2 <: Union{<:AbstractArray{RT, 2}, <:AbstractArray{RT, 3}},
#     R3 <: Union{<:AbstractArray{RT, 3}, <:AbstractArray{RT, 4}},
#     R4 <: Union{<:AbstractArray{RT, 4}, <:AbstractArray{RT, 5}}
# } <: AbstractH1OrL2Interpolants
#     ws::R1
#     ξs::R2
#     Ns::R2
#     ∇N_ξs::R3
#     ∇∇N_ξs::R4
# end 

# num_quadrature_points(interps::H1OrL2Interpolants) = size(interps.ws, 1)
# quadrature_point(interps::H1OrL2Interpolants, q::Int) = view(interps.ξs, :, q)
# quadrature_point(interps::H1OrL2Interpolants, q::Int, f::Int) = view(interps.ξs, :, q, f)
# quadrature_weight(interps::H1OrL2Interpolants, q::Int) = interps.ws[q]
# quadrature_weight(interps::H1OrL2Interpolants, q::Int, f::Int) = interps.ws[q, f]
# shape_function_gradient(interps::H1OrL2Interpolants, q::Int) = view(interps.∇N_ξs, :, :, q)
# shape_function_gradient(interps::H1OrL2Interpolants, q::Int, f::Int) = view(interps.∇N_ξs, :, :, q, f)
# shape_function_hessian(interps::H1OrL2Interpolants, q::Int) = view(interps.∇∇N_ξs, :, :, :, q)
# shape_function_hessian(interps::H1OrL2Interpolants, q::Int, f::Int) = view(interps.∇∇N_ξs, :, :, :, q, f)
# shape_function_value(interps::H1OrL2Interpolants, q::Int) = view(interps.Ns, :, q)
# shape_function_value(interps::H1OrL2Interpolants, q::Int, f::Int) = view(interps.Ns, :, q, f)

abstract type AbstractStaticH1OrL2Interpolants <: AbstractStaticInterpolants end
shape_function_gradient(interps::AbstractVector{T}, q::Int) where T <: AbstractStaticH1OrL2Interpolants = interps[q].∇N_ξ
shape_function_gradient(interps::AbstractMatrix{T}, q::Int, f::Int) where T <: AbstractStaticH1OrL2Interpolants = interps[q, f].∇N_ξ

struct StaticH1OrL2Interpolants{
    RT <: Number,
    ND, NN, NDNN
} <: AbstractStaticH1OrL2Interpolants
    w::RT
    ξ::SVector{ND, RT}
    N::SVector{NN, RT}
    ∇N_ξ::SMatrix{NN, ND, RT, NDNN}
end

function StaticH1OrL2Interpolants(el_type::AbstractElementType, X, ξ, w)
    ND = dimension(el_type)
    NN = num_cell_dofs(el_type)
    return StaticH1OrL2Interpolants(
        w, SVector{ND, Float64}(ξ),
        SVector{NN, Float64}(shape_function_value(el_type, X, ξ)),
        SMatrix{NN, ND, Float64, ND * NN}(shape_function_gradient(el_type, X, ξ))
    )
end

struct StaticH1OrL2InterpolantsWithHessians{
    RT <: Number,
    ND, NN, NDNN, NDNDNN
} <: AbstractStaticH1OrL2Interpolants
    w::RT
    ξ::SVector{ND, RT}
    N::SVector{NN, RT}
    ∇N_ξ::SMatrix{NN, ND, RT, NDNN}
    ∇∇N_ξ::SArray{Tuple{NN, ND, ND}, RT, 3, NDNDNN}
end

function StaticH1OrL2InterpolantsWithHessians(el_type::AbstractElementType, X, ξ, w)
    ND = dimension(el_type)
    NN = num_cell_dofs(el_type)
    return StaticH1OrL2InterpolantsWithHessians(
        w, SVector{ND, Float64}(ξ),
        SVector{NN, Float64}(shape_function_value(el_type, X, ξ)),
        SMatrix{NN, ND, Float64, ND * NN}(shape_function_gradient(el_type, X, ξ)),
        SArray{Tuple{NN, ND, ND}, Float64, 3, ND * ND * NN}(shape_function_hessian(el_type, X, ξ))
    )
end

shape_function_hessian(interps::AbstractVector{T}, q::Int) where T <: StaticH1OrL2InterpolantsWithHessians = interps[q].∇∇N_ξ
shape_function_hessian(interps::AbstractMatrix{T}, q::Int, f::Int) where T <: StaticH1OrL2InterpolantsWithHessians = interps[q, f].∇∇N_ξ

abstract type AbstractStaticHdivInterpolants <: AbstractStaticInterpolants end


struct StaticHdivInterpolants{
    RT <: Number,
    ND, NN, NDNN
} <: AbstractStaticHdivInterpolants
    w::RT
    ξ::SVector{ND, RT}
    N::SMatrix{NN, ND, RT, NDNN}
    divN_ξ::SVector{NN, RT}
    geometry_N::SVector{NN, RT}
    geometry_∇N_ξ::SMatrix{NN, ND, RT, NDNN}
end

function StaticHdivInterpolants(el_type::AbstractElementType, X, ξ, w)
    ND = dimension(el_type)
    NN = num_cell_dofs(el_type)
    return StaticHdivInterpolants(
        w, SVector{ND, eltype(ξ)}(ξ),
        SMatrix{NN, ND, eltype(ξ), ND * NN}(shape_function_value(el_type, X, ξ)),
        SVector{NN, eltype(ξ)}(shape_function_divergence(el_type, X, ξ)),
        SVector{NN, eltype(ξ)}(geometry_shape_function_value(el_type, X, ξ)),
        SMatrix{NN, ND, eltype(ξ), ND * NN}(geometry_shape_function_gradient(el_type, X, ξ))
    )
end

shape_function_divergence(interps::AbstractVector{T}, q::Int) where T <: StaticHdivInterpolants = interps[q].divN_ξ
shape_function_divergence(interps::AbstractVector{T}, q::Int, f::Int) where T <: StaticHdivInterpolants = interps[q, f].divN_ξ

# function _setup_cell_interpolants(
#     el_type::AbstractElementType{Lagrange, PD},
#     q_rule::AbstractQuadratureType,
#     ::Type{<:H1OrL2Interpolants}
# ) where PD
#     Xs = dof_coordinates(el_type)
#     ξs, ws = cell_quadrature_points_and_weights(el_type, q_rule)
    
#     ND = dimension(el_type)
#     NE = num_cell_dofs(el_type)
#     NQ = length(ws)

#     Ns = zeros(NE, NQ)
#     ∇N_ξs = zeros(ND, NE, NQ)
#     ∇∇N_ξs = zeros(ND, ND, NE, NQ)

#     for (q, ξ) in enumerate(eachcol(ξs))
#         if ND == 1
#             Ns[:, q] = shape_function_value(el_type, Xs, ξ[1])
#         else
#             Ns[:, q] = shape_function_value(el_type, Xs, ξ)
#         end
#         ∇N_ξs[:, :, q] = shape_function_gradient(el_type, Xs, ξ)
#         ∇∇N_ξs[:, :, :, q] = shape_function_hessian(el_type, Xs, ξ)
#     end

#     return H1OrL2Interpolants(ws, ξs, Ns, ∇N_ξs, ∇∇N_ξs)
# end

function _setup_cell_interpolants(
    el_type::AbstractElementType,
    q_rule::AbstractQuadratureType,
    type::Type{<:AbstractStaticInterpolants}
)
    Xs = dof_coordinates(el_type)
    ξs, ws = cell_quadrature_points_and_weights(el_type, q_rule)
    NN = num_cell_dofs(el_type)
    ND = dimension(el_type)
    if type <: StaticH1OrL2Interpolants
        interps = Vector{StaticH1OrL2Interpolants{Float64, ND, NN, ND * NN}}(undef, length(ws))
    elseif type <: StaticH1OrL2InterpolantsWithHessians
        interps = Vector{StaticH1OrL2InterpolantsWithHessians{Float64, ND, NN, ND * NN, ND * ND * NN}}(undef, length(ws))
    elseif type <: StaticHdivInterpolants
        interps = Vector{StaticHdivInterpolants{Float64, ND, NN, ND * NN}}(undef, length(ws))
    else
        @assert false "Unsupported type $type"
    end

    for (q, (w, ξ)) in enumerate(zip(ws, eachcol(ξs)))
        if typeof(el_type) <: Edge
            ξ = ξ[1]
        end
        interps[q] = type(el_type, Xs, ξ, w)
    end
    return interps
end

# function _setup_surface_interpolants(
#     el_type::AbstractElementType{Lagrange, PD},
#     q_rule::AbstractQuadratureType,
#     ::Type{<:H1OrL2Interpolants}
# ) where PD
#     if dimension(el_type) == 1
#         Xs = zeros(1, 1, 2)
#         Xs[1, 1, 1] = -1.
#         Xs[1, 1, 2] = 1.
#     else
#         Xs = dof_coordinates(el_type)[:, boundary_dofs(el_type)]
#     end
#     ξs, ws = surface_quadrature_points_and_weights(el_type, q_rule)
 
#     ND = dimension(el_type)
#     NE = num_cell_dofs(el_type)
#     NQ = size(ξs, 2)
#     NF = size(ξs, 3)

#     Ns = zeros(NE, NQ, NF)
#     ∇N_ξs = zeros(ND, NE, NQ, NF)
#     ∇∇N_ξs = zeros(ND, ND, NE, NQ, NF)

#     if dimension(el_type) == 1
#         Ns[1, :, 1] .= 1.
#         Ns[2, :, 2] .= 1.
#     else
#         for f in axes(ξs, 3)
#             for q in axes(ξs, 2)
#                 Ns[:, q, f] .= shape_function_value(el_type, Xs[:, :, f], ξs[:, q, f])
#                 ∇N_ξs[:, :, q, f] .= shape_function_gradient(el_type, Xs[:, :, f], ξs[:, q, f])
#                 ∇∇N_ξs[:, :, :, q, f] .= shape_function_hessian(el_type, Xs[:, :, f], ξs[:, q, f])
#             end
#         end
#     end
#     return H1OrL2Interpolants(ws, ξs, Ns, ∇N_ξs, ∇∇N_ξs)
# end

function _setup_surface_interpolants(
    el_type::AbstractElementType,
    q_rule::AbstractQuadratureType,
    type::Type{<:AbstractStaticInterpolants}
)
    if dimension(el_type) == 1
        Xs = zeros(1, 1, 2)
        Xs[1, 1, 1] = -1.
        Xs[1, 1, 2] = 1.
    else
        Xs = dof_coordinates(el_type)[:, boundary_dofs(el_type)]
    end
    ξs, ws = surface_quadrature_points_and_weights(el_type, q_rule)

    NN = num_cell_dofs(el_type)
    ND = dimension(el_type)
    NB = num_boundaries(el_type)
    if type <: StaticH1OrL2Interpolants
        interps = Matrix{type{Float64, ND, NN, ND * NN}}(undef, size(ws, 1), NB)
    elseif type <: StaticH1OrL2InterpolantsWithHessians
        interps = Matrix{type{Float64, ND, NN, ND * NN, ND * ND * NN}}(undef, size(ws, 1), NB)
    elseif type <: StaticHdivInterpolants
        interps = Matrix{type{Float64, ND, NN, ND * NN}}(undef, size(ws, 1), NB)
    else
        @assert false
    end

    # special annoying case that may be Lagrange specific
    if dimension(el_type) == 1
        if type <: StaticH1OrL2Interpolants
            interps[1, 1] = StaticH1OrL2Interpolants(
                1., SVector{ND, Float64}(ξs[1, 1, 1]),
                zero(SVector{NN, Float64}),
                zero(SMatrix{NN, ND, Float64, ND * NN})
            )
            interps[1, 2] = StaticH1OrL2Interpolants(
                1., SVector{ND, Float64}(ξs[1, 1, 2]),
                zero(SVector{NN, Float64}),
                zero(SMatrix{NN, ND, Float64, ND * NN})
            )
        elseif type <: StaticH1OrL2InterpolantsWithHessians
            interps[1, 1] = StaticH1OrL2InterpolantsWithHessians(
                1., SVector{ND, Float64}(ξs[1, 1, 1]),
                zero(SVector{NN, Float64}),
                zero(SMatrix{NN, ND, Float64, ND * NN}),
                zero(SArray{Tuple{NN, ND, ND}, Float64, 3, ND * ND * NN})
            )
            interps[1, 2] = StaticH1OrL2InterpolantsWithHessians(
                1., SVector{ND, Float64}(ξs[1, 1, 2]),
                zero(SVector{NN, Float64}),
                zero(SMatrix{NN, ND, Float64, ND * NN}),
                zero(SArray{Tuple{NN, ND, ND}, Float64, 3, ND * ND * NN})
            )
        else
            @assert false
        end
    else
        for f in axes(ξs, 3)
            for q in axes(ξs, 2)
                if typeof(el_type) <: Edge
                    ξ = ξs[1, q, f]
                else
                    ξ = ξs[:, q, f]
                end
                interps[q, f] = type(el_type, Xs[:, :, f], ξ, ws[q, f])
            end
        end
    end
    return interps
end

struct ReferenceFE{
    EType       <: AbstractElementType,
    BDofs       <: AbstractVector{<:SVector},
    BNorms      <: AbstractMatrix{<:Number},
    CellInterps <: Union{
        <:AbstractInterpolants, 
        <:AbstractVector{<:AbstractInterpolants}
    },
    SurfInterps <: Union{
        <:AbstractInterpolants,
        <:AbstractMatrix{<:AbstractInterpolants}
    }
}
    element::EType
    boundary_dofs::BDofs
    boundary_normals::BNorms
    cell_interps::CellInterps
    surf_interps::SurfInterps

    function ReferenceFE(
        el_type::AbstractElementType{Lagrange, PD},
        q_rule::AbstractQuadratureType;
        interpolants_type = StaticH1OrL2Interpolants
    ) where PD
        return ReferenceFE(el_type, q_rule, interpolants_type)
    end

    function ReferenceFE(
        el_type::AbstractElementType{RaviartThomas, PD},
        q_rule::AbstractQuadratureType;
        interpolants_type = StaticHdivInterpolants
    ) where PD
        return ReferenceFE(el_type, q_rule, interpolants_type)
    end

    function ReferenceFE(
        el_type::AbstractElementType,
        q_rule::AbstractQuadratureType,
        interpolants_type::Type{<:AbstractInterpolants}
    )
        bdofs = boundary_dofs(el_type)
        bdofs = map(x -> SVector{length(x), Int}(x), eachcol(bdofs))
        bnorms = boundary_normals(el_type)
        cell_interps = _setup_cell_interpolants(el_type, q_rule, interpolants_type)
        surf_interps = _setup_surface_interpolants(el_type, q_rule, interpolants_type)
        new{
            typeof(el_type), 
            typeof(bdofs), typeof(bnorms),
            typeof(cell_interps), typeof(surf_interps),
        }(
            el_type, bdofs, bnorms, cell_interps, surf_interps
        )
    end

    function ReferenceFE(
        element, bdofs, bnorms, cell_interps, surf_interps
    )
        new{
            typeof(element), 
            typeof(bdofs), typeof(bnorms),
            typeof(cell_interps), typeof(surf_interps)
        }(element, bdofs, bnorms, cell_interps, surf_interps)
    end
end

function Adapt.adapt_structure(to, re::ReferenceFE)
    bdofs = adapt(to, re.boundary_dofs)
    bnorms = adapt(to, re.boundary_normals)
    cell_interps = adapt(to, re.cell_interps)
    surf_interps = adapt(to, re.surf_interps)
    return ReferenceFE(re.element, bdofs, bnorms, cell_interps, surf_interps)
end

function Base.show(io::IO, re::ReferenceFE)
    println(io, "ReferenceFE")
    println(io, "  Element type      = $(element(re))")
    println(io, "  Polynomial type   = $(polynomial_type(element(re)))")
    println(io, "  Polynomial degree = $(polynomial_degree(element(re)))")
end

# topology interface
boundary_element(re::ReferenceFE, i::Int) = boundary_element(re.element, i)
boundary_normal(re::ReferenceFE, f::Int) = view(re.boundary_normals, :, f)
dimension(re::ReferenceFE) = dimension(re.element)
element(re::ReferenceFE) = re.element
num_boundaries(re::ReferenceFE) = num_boundaries(re.element)
num_edges(re::ReferenceFE) = num_edges(re.element)
num_faces(re::ReferenceFE) = num_faces(re.element)
polynomial_degree(re::ReferenceFE) = polynomial_degree(re.element)
vertex_coordinates(re::ReferenceFE) = vertex_coordinates(re.element)

# dof interface
boundary_dofs(re::ReferenceFE) = re.boundary_dofs
boundary_dofs(re::ReferenceFE, f::Int) = re.boundary_dofs[f]
dof_coordinates(re::ReferenceFE) = dof_coordinates(re.element)
interior_dofs(re::ReferenceFE) = interior_dofs(re.element)
num_cell_dofs(re::ReferenceFE) = num_cell_dofs(re.element)

# quadrature interface
cell_quadrature_point(re::ReferenceFE, q::Int) = quadrature_point(re.cell_interps, q)
cell_quadrature_weight(re::ReferenceFE, q::Int) = quadrature_weight(re.cell_interps, q)
num_cell_quadrature_points(re::ReferenceFE) = num_quadrature_points(re.cell_interps)
num_surface_quadrature_points(re::ReferenceFE) = num_quadrature_points(re.surf_interps)
surface_quadrature_point(re::ReferenceFE, q::Int, f::Int) = quadrature_point(re.surf_interps, q, f)
surface_quadrature_weight(re::ReferenceFE, q::Int, f::Int) = quadrature_weight(re.surf_interps, q, f)

# interpolation interface
cell_shape_function_divergence(re::ReferenceFE, q::Int) = shape_function_divergence(re.cell_interps, q)
surface_shape_function_divergence(re::ReferenceFE, q::Int, f::Int) = shape_function_divergence(re.surf_interps, q, f)
cell_shape_function_gradient(re::ReferenceFE, q::Int) = shape_function_gradient(re.cell_interps, q)
surface_shape_function_gradient(re::ReferenceFE, q::Int, f::Int) = shape_function_gradient(re.surf_interps, q, f)
cell_shape_function_hessian(re::ReferenceFE, q::Int) = shape_function_hessian(re.cell_interps, q)
surface_shape_function_hessian(re::ReferenceFE, q::Int, f::Int) = shape_function_hessian(re.surf_interps, q, f)
cell_shape_function_value(re::ReferenceFE, q::Int) = shape_function_value(re.cell_interps, q)
surface_shape_function_value(re::ReferenceFE, q::Int, f::Int) = shape_function_value(re.surf_interps, q, f)

struct MappedH1OrL2Interpolants{A, B, C, D}
    X_q::A
    N::B
    ∇N_X::C
    JxW::D
end

function MappedH1OrL2Interpolants(e::ReferenceFE, X, q)
    # unpacks
    w = cell_quadrature_weight(e, q)
    N = cell_shape_function_value(e, q)
    ∇N_ξ = cell_shape_function_gradient(e, q)
  
    # for N x D x D behavior
    # interpolate coordinates
    X_q = X * N

    # map shape function gradients
    J = (X * ∇N_ξ)'
    J_inv = inv(J)
    ∇N_X = (J_inv * ∇N_ξ')'

    # JxW
    JxW = det(J) * w
  
    return MappedH1OrL2Interpolants(X_q, N, ∇N_X, JxW)
end

function MappedH1OrL2Interpolants(interps::StaticH1OrL2Interpolants, X)
    w, N, ∇N_ξ = interps.w, interps.N, interps.∇N_ξ
  
    # interpolate coordinates
    X_q = X * N
    
    # N x D x D behavior
    # map shape function gradients
    J = (X * ∇N_ξ)'
    J_inv = inv(J)
    ∇N_X = (J_inv * ∇N_ξ')'

    # JxW
    JxW = det(J) * w
  
    return MappedH1OrL2Interpolants(X_q, N, ∇N_X, JxW)
end

dimension(interps::MappedH1OrL2Interpolants) = size(interps.∇N_X, 1)

struct MappedH1OrL2SurfaceInterpolants{
    A <: AbstractArray, 
    B <: AbstractArray, 
    C <: AbstractArray,
    D <: Number,
    E <: AbstractArray
  }
    X_q::A
    N::B
    N_reduced::C
    JxW::D
    n::E
end

# specialize for surface shape functions
function MappedH1OrL2SurfaceInterpolants(e::ReferenceFE, X, q::Integer, f::Integer)
    w = surface_quadrature_weight(e, q, f)
    N = surface_shape_function_value(e, q, f)
    NNPS = num_cell_dofs(boundary_element(e, f))
    edge_nodes = SVector{NNPS, Int}(boundary_dofs(e, f))
    N_reduced = SVector{NNPS, eltype(N)}(@views N[edge_nodes])
    n = boundary_normal(e, f)

    # jacobian
    X_diff = X[:, 2] - X[:, 1]
    det_J = norm(X_diff)
    # interpolate coordinates
    edge_nodes = boundary_dofs(e, f)
    X_q = SVector{2, eltype(X)}(@views X[:, edge_nodes] * N_reduced)
  
    # JxW
    JxW = det_J * w
  
    # TODO below incorrect. Not giving correct gradient
    # or normal
    # @show N
    return MappedH1OrL2SurfaceInterpolants(X_q, N, N_reduced, JxW, n)
end

@inline function _piola_map(
    J::SMatrix{2, 2, T, 4},
    detJ::T,
    N::SMatrix{ND, 2, T, ND2}
) where {T, ND, ND2}

    invdetJ = inv(detJ)

    return SMatrix{ND, 2, T, ND2}(
        ntuple(i -> begin
            a = (i - 1) ÷ 2 + 1   # row
            j = (i - 1) % 2 + 1   # component
            (J * N[a, :])[j] * invdetJ
        end, ND2)
    )
end

struct MappedHdivInterpolants{A, B, C, D}
    X_q::A
    N::B
    divN_X::C
    JxW::D
end

function MappedHdivInterpolants(
    interps::StaticHdivInterpolants, X, orientations
)
    w, N_ξ, divN_ξ = interps.w, interps.N, interps.divN_ξ
    N = interps.geometry_N
    ∇N_ξ = interps.geometry_∇N_ξ

    # interpolate coordinates
    X_q = X * N

    # calculate J using geometry shape function grad (H1)
    J = (X * ∇N_ξ)'
    detJ = det(J)

    # map shape function values using piola transform
    N = _piola_map(J, detJ, N_ξ)
    N = orientations .* N

    # map divergence
    divN_X = divN_ξ / detJ

    # JxW
    JxW = detJ * w
    return MappedHdivInterpolants(X_q, N, divN_X, JxW)
end
