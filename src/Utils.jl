function types_to_generate_quadrature(e::E) where E <: ReferenceFEType
  n_nodes, dim = num_nodes(e), num_dimensions(e)

  eltypes = [Float32, Float64]

  types   = ()
  for n in axes(eltypes, 1)
    type  = eltypes[n]
    types = (types...,
      (
        :($type),
        :(SVector{$dim, $type}),
      )
    )
    types = (types...,
    (
      :($type),
      :(MVector{$dim, $type}),
    )
    )
  end
  return types
end

function types_to_generate_interpolants(e::E) where E <: ReferenceFEType

  types = ()
  types = (types...,
    (
      :(SVector),
      :(SMatrix),
      :(SArray)
    )
  )
  types = (types...,
    (
      :(MVector),
      :(MMatrix),
      :(MArray)
    )
  )
  return types
end