struct H1OrL2CellInterpolants{
    T <: Number,
    W <: AbstractVector{T},
    V <: AbstractMatrix{T},
    G <: AbstractArray{T, 3},
    H <: AbstractArray{T, 4}
} <: AbstractDynamicInterpolants{T}
    weights::W
    values::V
    gradients::G
    hessians::H
end

function Adapt.adapt_structure(to, interps::H1OrL2CellInterpolants)
    return H1OrL2CellInterpolants(
        adapt(to, interps.weights),
        adapt(to, interps.values),
        adapt(to, interps.gradients),
        adapt(to, interps.hessians)
    )
end

struct H1OrL2SurfaceInterpolants{
    T <: Number,
    W <: AbstractMatrix{T},
    V <: AbstractArray{T, 3},
    G <: AbstractArray{T, 4},
    H <: AbstractArray{T, 5}
} <: AbstractDynamicInterpolants{T}
    weights::W
    values::V
    gradients::G
    hessians::H
end

function Adapt.adapt_structure(to, interps::H1OrL2SurfaceInterpolants)
    return H1OrL2SurfaceInterpolants(
        adapt(to, interps.weights),
        adapt(to, interps.values),
        adapt(to, interps.gradients),
        adapt(to, interps.hessians)
    )
end

function _setup_cell_interpolants(
    el_type::AbstractElementType,
    q_rule::AbstractQuadratureType,
    ::Type{<:H1OrL2CellInterpolants},
    ::Type{T} = Float64
) where T
    NN = num_cell_dofs(el_type)
    ND = dimension(el_type)
    ξs, ws = cell_quadrature_points_and_weights(el_type, q_rule)
    Ns = Matrix{T}(undef, NN, length(ws))
    ∇N_ξs = Array{T, 3}(undef, ND, NN, length(ws))
    ∇∇N_ξs = Array{T, 4}(undef, ND, ND, NN, length(ws))
    for q in axes(ws, 1)
        ξ = SVector{dimension(el_type), T}(@views ξs[:, q])
        Ns[:, q] = shape_function_value(el_type, ξ)
        ∇N_ξs[:, :, q] = shape_function_gradient(el_type, ξ)'
        ∇∇N_temp = shape_function_hessian(el_type, ξ)
        for n in axes(∇∇N_temp, 1)
            for d1 in axes(∇∇N_temp, 2)
                for d2 in axes(∇∇N_temp, 3)
                    ∇∇N_ξs[d2, d1, n, q] = ∇∇N_temp[n, d1, d2]
                end
            end
        end
    end
    return H1OrL2CellInterpolants(ws, Ns, ∇N_ξs, ∇∇N_ξs)
end

function _setup_surface_interpolants(
    el_type::AbstractElementType,
    q_rule::AbstractQuadratureType,
    ::Type{<:H1OrL2SurfaceInterpolants},
    ::Type{T} = Float64
) where T <: Number
    NN = num_cell_dofs(el_type)
    ND = dimension(el_type)
    ξs, ws = surface_quadrature_points_and_weights(el_type, q_rule)
    Ns = Array{T, 3}(undef, num_cell_dofs(el_type), size(ws, 1), size(ws, 2))
    ∇N_ξs = Array{T, 4}(undef, ND, NN, size(ws, 1), size(ws, 2))
    ∇∇N_ξs = Array{T, 5}(undef, ND, ND, NN, size(ws, 1), size(ws, 2))
    for f in axes(ws, 2)
        for q in axes(ws, 1)
            ξ = SVector{dimension(el_type), T}(@views ξs[:, q, f])
            Ns[:, q, f] = shape_function_value(el_type, ξ)
            ∇N_ξs[:, :, q, f] = shape_function_gradient(el_type, ξ)'
            ∇∇N_temp = shape_function_hessian(el_type, ξ)
            for n in axes(∇∇N_temp, 1)
                for d1 in axes(∇∇N_temp, 2)
                    for d2 in axes(∇∇N_temp, 3)
                        ∇∇N_ξs[d2, d1, n, q, f] = ∇∇N_temp[n, d1, d2]
                    end
                end
            end
        end
    end
    return H1OrL2SurfaceInterpolants(ws, Ns, ∇N_ξs, ∇∇N_ξs)
end

# TODO add some info about e.g. element type
# but do it in a way that is type stable for all
# blocks in a function space
struct ReferenceFE_v2{
    D,
    I <: Integer,
    T <: Number,
    B <: AbstractMatrix{I},
    C <: Union{<:AbstractInterpolants{T}, <:AbstractVector{<:AbstractInterpolants{T}}},
    S <: Union{<:AbstractInterpolants{T}, <:AbstractMatrix{<:AbstractInterpolants{T}}},
}
    boundary_dofs::B
    boundary_normals::B
    cell_interps::C
    surf_interps::S

    function ReferenceFE_v2{D, I, T}(bdofs::B, bnorms, cinterps::C, sinterps::S) where {D, I, T, B <: AbstractMatrix{I}, C, S}
        new{D, I, T, B, C, S}(bdofs, bnorms, cinterps, sinterps)
    end
end

_cell_interpolants_type(
    ::AbstractElementType{D, Lagrange}, ::Type{<:Array}
) where D = H1OrL2CellInterpolants
_cell_interpolants_type(
    ::AbstractElementType{D, Lagrange}, ::Type{<:SArray}
) where D = StaticH1OrL2Interpolants
_surface_interpolants_type(
    ::AbstractElementType{D, Lagrange}, ::Type{<:Array}
) where D = H1OrL2SurfaceInterpolants
_surface_interpolants_type(
    ::AbstractElementType{D, Lagrange}, ::Type{<:SArray}
) where D = StaticH1OrL2Interpolants

# need to specialize on space type
function ReferenceFE_v2(
    el_type::AbstractElementType,
    q_rule::AbstractQuadratureType,
    float_type::Type{T} = Float64;
    array_type::Type = SArray
) where T <: Number
    bdofs = boundary_dofs(el_type)
    bnorms = boundary_normals(el_type)
    cell_interps = _setup_cell_interpolants(el_type, q_rule, _cell_interpolants_type(el_type, array_type), float_type)
    surf_interps = _setup_surface_interpolants(el_type, q_rule, _surface_interpolants_type(el_type, array_type), float_type)
    return ReferenceFE_v2{dimension(el_type), Int, T}(bdofs, bnorms, cell_interps, surf_interps)
end

function Adapt.adapt_structure(to, ref_fe::ReferenceFE_v2{D, I, T}) where {D, I, T}
    return ReferenceFE_v2{D, I, T}(
        adapt(to, ref_fe.boundary_dofs),
        adapt(to, ref_fe.boundary_normals),
        adapt(to, ref_fe.cell_interps),
        adapt(to, ref_fe.surf_interps)
    )
end
