struct Vertex <: AbstractVertex
end

function shape_function_value(::Vertex, _, _)
    return [1.]
end

function shape_function_gradient(::Vertex, _, _)
    return Matrix{Float64}(undef, 0, 1)
end

function shape_function_hessian(::Vertex, _, _)
    return Array{Float64, 3}(undef, 0, 0, 1)
end
