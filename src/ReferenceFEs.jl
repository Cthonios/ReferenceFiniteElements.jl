abstract type AbstractInterpolants end

struct H1OrL2Interpolants{
    RT <: Number,
    R1 <: Union{<:AbstractArray{RT, 1}, <:AbstractArray{RT, 2}},
    R2 <: Union{<:AbstractArray{RT, 2}, <:AbstractArray{RT, 3}},
    R3 <: Union{<:AbstractArray{RT, 3}, <:AbstractArray{RT, 4}},
    R4 <: Union{<:AbstractArray{RT, 4}, <:AbstractArray{RT, 5}}
} <: AbstractInterpolants
    ws::R1
    ξs::R2
    Ns::R2
    ∇N_ξs::R3
    ∇∇N_ξs::R4
end 

num_quadrature_points(interps::H1OrL2Interpolants) = size(interps.ws, 1)
quadrature_point(interps::H1OrL2Interpolants, q::Int) = view(interps.ξs, :, q)
quadrature_point(interps::H1OrL2Interpolants, q::Int, f::Int) = view(interps.ξs, :, q, f)
quadrature_weight(interps::H1OrL2Interpolants, q::Int) = interps.ws[q]
quadrature_weight(interps::H1OrL2Interpolants, q::Int, f::Int) = interps.ws[q, f]
shape_function_gradient(interps::H1OrL2Interpolants, q::Int) = view(interps.∇N_ξs, :, :, q)
shape_function_gradient(interps::H1OrL2Interpolants, q::Int, f::Int) = view(interps.∇N_ξs, :, :, q, f)
shape_function_hessian(interps::H1OrL2Interpolants, q::Int) = view(interps.∇∇N_ξs, :, :, :, q)
shape_function_hessian(interps::H1OrL2Interpolants, q::Int, f::Int) = view(interps.∇∇N_ξs, :, :, :, q, f)
shape_function_value(interps::H1OrL2Interpolants, q::Int) = view(interps.Ns, :, q)
shape_function_value(interps::H1OrL2Interpolants, q::Int, f::Int) = view(interps.Ns, :, q, f)


struct SH1OrL2Interpolants{
    RT <: Number,
    ND, NN, NDNN, NDNDNN
} <: AbstractInterpolants
    w::RT
    ξ::SVector{ND, RT}
    N::SVector{NN, RT}
    ∇N_ξ::SMatrix{ND, NN, RT, NDNN}
    ∇∇N_ξ::SArray{Tuple{ND, ND, NN}, RT, 3, NDNDNN}
end

num_quadrature_points(interps::AbstractVector{T}) where T <: SH1OrL2Interpolants = length(interps)
num_quadrature_points(interps::AbstractMatrix{T}) where T <: SH1OrL2Interpolants = size(interps, 1)
quadrature_point(interps::AbstractVector{T}, q::Int) where T <: SH1OrL2Interpolants = interps[q].ξ
quadrature_point(interps::AbstractMatrix{T}, q::Int, f::Int) where T <: SH1OrL2Interpolants = interps[q, f].ξ
quadrature_weight(interps::AbstractVector{T}, q::Int) where T <: SH1OrL2Interpolants = interps[q].w
quadrature_weight(interps::AbstractMatrix{T}, q::Int, f::Int) where T <: SH1OrL2Interpolants = interps[q, f].w

function _setup_cell_interpolants(
    el_type::AbstractElementType{Lagrange, PD},
    q_rule::AbstractQuadratureType,
    ::Type{<:H1OrL2Interpolants}
) where PD
    Xs = dof_coordinates(el_type)
    ξs, ws = cell_quadrature_points_and_weights(el_type, q_rule)
    
    ND = dimension(el_type)
    NE = num_cell_dofs(el_type)
    NQ = length(ws)

    Ns = zeros(NE, NQ)
    ∇N_ξs = zeros(ND, NE, NQ)
    ∇∇N_ξs = zeros(ND, ND, NE, NQ)

    for (q, ξ) in enumerate(eachcol(ξs))
        if ND == 1
            Ns[:, q] = shape_function_value(el_type, Xs, ξ[1])
        else
            Ns[:, q] = shape_function_value(el_type, Xs, ξ)
        end
        ∇N_ξs[:, :, q] = shape_function_gradient(el_type, Xs, ξ)
        ∇∇N_ξs[:, :, :, q] = shape_function_hessian(el_type, Xs, ξ)
    end

    return H1OrL2Interpolants(ws, ξs, Ns, ∇N_ξs, ∇∇N_ξs)
end

function _setup_cell_interpolants(
    el_type::AbstractElementType{Lagrange, PD},
    q_rule::AbstractQuadratureType,
    ::Type{<:SH1OrL2Interpolants}
) where PD
    Xs = dof_coordinates(el_type)
    ξs, ws = quadrature_points_and_weights(el_type, q_rule)
    NN = num_cell_dofs(el_type)
    ND = dimension(el_type)
    interps = Vector{SH1OrL2Interpolants{Float64, ND, NN, ND * NN, ND * ND * NN}}(undef, length(ws))
    map!(
        (w, ξ) -> SH1OrL2Interpolants(
            w, 
            SVector{ND, Float64}(ξ),
            SVector{NN, Float64}(shape_function_value(el_type, Xs, ξ)),
            SMatrix{ND, NN, Float64, ND * NN}(shape_function_gradient(el_type, Xs, ξ)),
            SArray{Tuple{ND, ND, NN}, Float64, 3, ND * ND * NN}(shape_function_hessian(el_type, Xs, ξ)) 
        ),
        interps,
        ws, eachcol(ξs)
    )
    return interps
end

function _setup_surface_interpolants(
    el_type::AbstractElementType{Lagrange, PD},
    q_rule::AbstractQuadratureType,
    ::Type{<:H1OrL2Interpolants}
) where PD
    if dimension(el_type) == 1
        Xs = zeros(1, 1, 2)
        Xs[1, 1, 1] = -1.
        Xs[1, 1, 2] = 1.
    else
        Xs = dof_coordinates(el_type)[:, edge_vertices(el_type)]
    end
    ξs, ws = surface_quadrature_points_and_weights(el_type, q_rule)
 
    ND = dimension(el_type)
    NE = num_cell_dofs(el_type)
    NQ = size(ξs, 2)
    NF = size(ξs, 3)

    Ns = zeros(NE, NQ, NF)
    ∇N_ξs = zeros(ND, NE, NQ, NF)
    ∇∇N_ξs = zeros(ND, ND, NE, NQ, NF)

    if dimension(el_type) == 1
        Ns[1, :, 1] .= 1.
        Ns[2, :, 2] .= 1.
    else
        for f in axes(ξs, 3)
            for q in axes(ξs, 2)
                Ns[:, q, f] .= shape_function_value(el_type, Xs[:, :, f], ξs[:, q, f])
                ∇N_ξs[:, :, q, f] .= shape_function_gradient(el_type, Xs[:, :, f], ξs[:, q, f])
                ∇∇N_ξs[:, :, :, q, f] .= shape_function_hessian(el_type, Xs[:, :, f], ξs[:, q, f])
            end
        end
    end
    return H1OrL2Interpolants(ws, ξs, Ns, ∇N_ξs, ∇∇N_ξs)
end

struct ReferenceFE{
    EType <: AbstractElementType,
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
    cell_interps::CellInterps
    surf_interps::SurfInterps

    function ReferenceFE(
        el_type::AbstractElementType,
        q_rule::AbstractQuadratureType;
        type = H1OrL2Interpolants
    )
        cell_interps = _setup_cell_interpolants(el_type, q_rule, type)
        surf_interps = _setup_surface_interpolants(el_type, q_rule, type)
        new{typeof(el_type), typeof(cell_interps), typeof(surf_interps)}(
            el_type, cell_interps, surf_interps
        )
    end
end

# topology interface
boundary_element(re::ReferenceFE) = boundary_element(re.element)
dimension(re::ReferenceFE) = dimension(re.element)
num_boundaries(re::ReferenceFE) = num_boundaries(re.element)
num_edges(re::ReferenceFE) = num_edges(re.element)
num_faces(re::ReferenceFE) = num_faces(re.element)
# num_facets(re::ReferenceFE) = num_facets(re.element)
polynomial_degree(re::ReferenceFE) = polynomial_degree(re.element)
# surface_element(re::ReferenceFE) = surface_element(re.element)
vertex_coordinates(re::ReferenceFE) = vertex_coordinates(re.element)

# dof interface
dof_coordinates(re::ReferenceFE) = dof_coordinates(re.element)
interior_dofs(re::ReferenceFE) = interior_dofs(re.element)

# quadrature interface
cell_quadrature_point(re::ReferenceFE, q::Int) = quadrature_point(re.cell_interps, q)
cell_quadrature_weight(re::ReferenceFE, q::Int) = quadrature_weight(re.cell_interps, q)
num_cell_quadrature_points(re::ReferenceFE) = num_quadrature_points(re.cell_interps)
num_surface_quadrature_points(re::ReferenceFE) = num_quadrature_points(re.surf_interps)
surface_quadrature_point(re::ReferenceFE, q::Int, f::Int) = quadrature_point(re.surf_interps, q, f)
surface_quadrature_weight(re::ReferenceFE, q::Int, f::Int) = quadrature_weight(re.surf_interps, q, f)

# interpolation interface
cell_shape_function_gradient(re::ReferenceFE, q::Int) = shape_function_gradient(re.cell_interps, q)
surface_shape_function_gradient(re::ReferenceFE, q::Int, f::Int) = shape_function_gradient(re.surf_interps, q, f)
cell_shape_function_hessian(re::ReferenceFE, q::Int) = shape_function_hessian(re.cell_interps, q)
surface_shape_function_hessian(re::ReferenceFE, q::Int, f::Int) = shape_function_hessian(re.surf_interps, q, f)
cell_shape_function_value(re::ReferenceFE, q::Int) = shape_function_value(re.cell_interps, q)
surface_shape_function_value(re::ReferenceFE, q::Int, f::Int) = shape_function_value(re.surf_interps, q, f)