struct ArrayBackend{T <: AbstractArray}
end

get_backend(::ArrayBackend{T}) where T = T

# this one below is still floating around
convert_to(a::A, ::ArrayBackend{Array}) where A <: Real = [a]
convert_to(a::A, ::ArrayBackend{MArray}) where A <: Real = @MArray [a]
convert_to(a::A, ::ArrayBackend{SArray}) where A <: Real = @SArray [a]

convert_to(a::A, ::ArrayBackend{Array}) where A <: Vector = a
convert_to(a::A, ::ArrayBackend{MArray}) where A <: Vector = @MVector [a]
convert_to(a::A, ::ArrayBackend{SArray}) where A <: Vector = @SVector [a]

convert_to(a::A, ::ArrayBackend{Array}) where A <: Matrix = a
convert_to(a::A, ::ArrayBackend{MArray}) where A <: Matrix = @MMatrix [a]
convert_to(a::A, ::ArrayBackend{SArray}) where A <: Matrix = @SMatrix [a]

convert_to(::ArrayBackend{Array}, a...) = [a...]
convert_to(::ArrayBackend{MArray}, a...) = MVector(a)
convert_to(::ArrayBackend{SArray}, a...) = SVector(a)

function convert_to_3d_array(e::AbstractElementType{I, D, P, Q}, ::ArrayBackend{Array}, vals...) where {I, D, P, Q}
  ret = Array{eltype(vals), 3}(undef, num_shape_functions(e), D, D)
  n = 1
  for k in axes(ret, 3)
    for j in axes(ret, 2)
      for i in axes(ret, 1)
        ret[i, j, k] = vals[n]
        n = n + 1
      end
    end
  end
  return ret
end

function convert_to_3d_array(e::AbstractElementType{I, D, P, Q}, ::ArrayBackend{MArray}, vals...) where {I, D, P, Q}
  N = num_shape_functions(e)
  return MArray{Tuple{N, D, D}, eltype(vals), 3, N * D * D}(vals...)
end

function convert_to_3d_array(e::AbstractElementType{I, D, P, Q}, ::ArrayBackend{SArray}, vals...) where {I, D, P, Q}
  N = num_shape_functions(e)
  return SArray{Tuple{N, D, D}, eltype(vals), 3, N * D * D}(vals...)
end

function convert_to_matrix(e::AbstractElementType{I, D, P, Q}, ::ArrayBackend{Array}, vals...) where {I, D, P, Q}
  ret = Matrix{eltype(vals)}(undef, num_shape_functions(e), D)
  n = 1
  for j in axes(ret, 2)
    for i in axes(ret, 1)
      ret[i, j] = vals[n]
      n = n + 1
    end
  end
  return ret
  # return eltype(vals)[vals...]
end

function convert_to_matrix(e::AbstractElementType{I, D, P, Q}, ::ArrayBackend{MArray}, vals...) where {I, D, P, Q}
  N = num_shape_functions(e)
  return MMatrix{N, D, eltype(vals), N * D}(vals...)
end

function convert_to_matrix(e::AbstractElementType{I, D, P, Q}, ::ArrayBackend{SArray}, vals...) where {I, D, P, Q}
  N = num_shape_functions(e)
  return SMatrix{N, D, eltype(vals), N * D}(vals...)
end

function convert_to_vector(::AbstractElementType{I, D, P, Q}, ::ArrayBackend{Array}, vals...) where {I, D, P, Q}
  return eltype(vals)[vals...]
end

function convert_to_vector(e::AbstractElementType{I, D, P, Q}, ::ArrayBackend{MArray}, vals...) where {I, D, P, Q}
  N = num_shape_functions(e)
  return MVector{N, eltype(vals)}(vals...)
end

function convert_to_vector(e::AbstractElementType{I, D, P, Q}, ::ArrayBackend{SArray}, vals...) where {I, D, P, Q}
  N = num_shape_functions(e)
  return SVector{N, eltype(vals)}(vals...)
end

function convert_to_vector_coords(::AbstractElementType{I, D, P, Q}, ::ArrayBackend{Array}, vals...) where {I, D, P, Q}
  return eltype(vals)[vals...]
end

function convert_to_vector_coords(::AbstractElementType{I, D, P, Q}, ::ArrayBackend{MArray}, vals...) where {I, D, P, Q}
  return MVector{D, eltype(vals)}(vals...)
end

function convert_to_vector_coords(::AbstractElementType{I, D, P, Q}, ::ArrayBackend{SArray}, vals...) where {I, D, P, Q}
  return SVector{D, eltype(vals)}(vals...)
end
