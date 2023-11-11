module ReferenceFiniteElementsExodusExt

using Exodus
using ReferenceFiniteElements

name_to_type = Dict{String, Type{<:ReferenceFiniteElements.ReferenceFEType}}(
  "HEX8"    => Hex8,
  "TETRA4"  => Tet4,
  "TETRA10" => Tet10,
  "TRI3"    => Tri3,
  "TRI6"    => Tri6,
  "QUAD4"   => Quad4,
  "QUAD9"   => Quad9
)

function ReferenceFiniteElements.ReferenceFE(block::B, q_order::Int) where B <: Exodus.Block
  return ReferenceFiniteElements.ReferenceFE(name_to_type[block.elem_type](q_order))
end

end # module