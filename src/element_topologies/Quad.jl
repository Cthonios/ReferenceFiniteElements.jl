

function edges(e::AbstractQuad{V}, backend::ArrayBackend) where V
  edge = surface_element(e)
  edges = [
    convert_to(backend, 1, 5:nvertices(edge) - 2, 2)
    
  ]
end

function faces(::AbstractQuad{V}, backend::ArrayBackend) where V
  [convert_to(backend, 1:V)]
end

function vertices(::AbstractQuad{V}, backend::ArrayBackend) where V
  convert_to(backend, 1:V...)
end

function vertex_coordinates(e::AbstractQuad{V}, backend::ArrayBackend; exclude_interior_vertices=false) where V
  edge = surface_element(e)
  edge_coords = vertex_coordinates(edge, backend)
  coords = [
    convert_to_vector_coords(e, backend, edge_coords[1][1], -1.),
    convert_to_vector_coords(e, backend, 1., edge_coords[1][1]),
    convert_to_vector_coords(e, backend, edge_coords[end][1], 1.),
    convert_to_vector_coords(e, backend, -1., edge_coords[end][1]),
  ]

  if V > 4
    for val in edge_coords[2:end - 1]
      push!(coords, convert_to_vector_coords(e, backend, val[1], -1.))
    end

    for val in edge_coords[2:end - 1]
      push!(coords, convert_to_vector_coords(e, backend, 1., val[1]))
    end

    for val in edge_coords[2:end - 1]
      push!(coords, convert_to_vector_coords(e, backend, val[1], 1.))
    end

    for val in edge_coords[2:end - 1]
      push!(coords, convert_to_vector_coords(e, backend, -1., val[1]))
    end

    if !exclude_interior_vertices
      for val in Iterators.product(edge_coords[2:end - 1], edge_coords[2:end - 1])
        push!(coords, convert_to_vector_coords(e, backend, val[1][1], val[2][1]))
      end
    end
  end

  return coords
end

function normals(e::AbstractQuad{V}, backend::ArrayBackend) where V
  return [
    convert_to_vector_coords(e, backend, 0., -1.),
    convert_to_vector_coords(e, backend, 1., 0.),
    convert_to_vector_coords(e, backend, 0., 1.),
    convert_to_vector_coords(e, backend, -1., 0.),
  ]
end

# Need error checks for supported orders

struct Quad{V} <: AbstractQuad{V}
  function Quad{V}(; exclude_interior_vertices=false) where V
    check_1 = V isa Integer
    if exclude_interior_vertices
      @assert false "Not implemented yet"
    else
      check_2 = false
      for n in 2:V
        check_2 = check_2 || V == n * n
      end
      @assert check_1 && check_2 "Supported parameters are e.g. 4, 9, 16, etc."
    end
    # TODO more error checks for not including interiors
    new{V}()
  end
end

function surface_element(::Quad{V}) where V
  V_edge = isqrt(V)
  if V_edge * V_edge == V
    return Edge{V_edge}()
  else
    @assert "Not implemented yet"
  end
end
