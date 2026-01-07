module ExodusExt

using Exodus
using ReferenceFiniteElements

const name_to_type = Dict{String, Any}(
  "HEX"     => Hex,
  "HEX8"    => Hex,
  "TETRA"   => Tet,
  "TETRA4"  => Tet,
  "TETRA10" => Tet,
  "TRI"     => Tet,
  "TRI3"    => Tri,
  "TRI6"    => Tri,
  "QUAD"    => Quad,
  "QUAD4"   => Quad,
  "QUAD9"   => Quad
)
const name_to_p = Dict{String, Int}(
  "HEX"     => 1,
  "HEX8"    => 1,
  "TETRA"   => 1,
  "TETRA4"  => 1,
  "TETRA10" => 2,
  "TRI"     => 1,
  "TRI3"    => 1,
  "TRI6"    => 2,
  "QUAD"    => 1,
  "QUAD4"   => 1,
  "QUAD9"   => 2
)
const default_q = Dict{String, Int}(
  "HEX"     => 2,
  "HEX8"    => 2,
  "TETRA"   => 2,
  "TETRA4"  => 2,
  "TETRA10" => 2,
  "TRI"     => 2,
  "TRI3"    => 2,
  "TRI6"    => 2,
  "QUAD"    => 2,
  "QUAD4"   => 2,
  "QUAD9"   => 2
)

function ReferenceFiniteElements.ReferenceFE(
  block::B, 
  interp_type::Type{<:ReferenceFiniteElements.AbstractPolynomialType};
  p_order::Union{Int, Nothing} = nothing, 
  q_order::Union{Int, Nothing} = nothing
) where B <: Exodus.Block
  el_type = name_to_type[block.elem_type]

  if p_order === nothing
    p_order = name_to_p[block.elem_type]
  end

  if q_order === nothing
    q_order = default_q[block.elem_type]
  end

  return ReferenceFiniteElements.ReferenceFE(
    el_type{interp_type, p_order}(),
    GaussLobattoLegendre(q_order)
  )
end

end # module
