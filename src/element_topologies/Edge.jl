struct Edge{V} <: AbstractEdge{V}
  function Edge{V}() where V
    @assert V isa Integer && V > 1
    new{V}()
  end
end

function edges(::Edge{V}, backend::ArrayBackend) where V
  [convert_to(backend, 1:V...)]
end

function faces(::Edge{V}, backend::ArrayBackend) where V
  []  
end

function vertices(::Edge{V}, backend::ArrayBackend) where V 
  convert_to(backend, 1:V...)
end

function vertex_coordinates(e::Edge{V}, backend::ArrayBackend; interval=(-1., 1.)) where V
  return map(x -> convert_to_vector_coords(e, backend, x), LinRange(interval[1], interval[2], V))
end

function normals(e::Edge, backend::ArrayBackend) 
  return map(x -> convert_to_vector_coords(e, backend, x), LinRange(-1., 1., 2))
end
