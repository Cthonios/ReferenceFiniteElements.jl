module ReferenceFiniteElementsExodusExt

using Exodus
using ReferenceFiniteElements
using StaticArrays

name_to_type = Dict{String, Type{<:ReferenceFiniteElements.ReferenceFEType}}(
  "HEX8"    => Hex8,
  "TETRA4"  => Tet4,
  "TETRA10" => Tet10,
  "TRI3"    => Tri3,
  "TRI6"    => Tri6,
  "QUAD4"   => Quad4,
  "QUAD9"   => Quad9
)

function ReferenceFiniteElements.ReferenceFE(
  block::B, q_order::Int;
  int_type   = Int64,
  float_type = Float64,
  array_type = SArray
) where B <: Exodus.Block
  return ReferenceFiniteElements.ReferenceFE(
    name_to_type[block.elem_type](Val(q_order));
    int_type=int_type, float_type=float_type,
    array_type=array_type
  )
end

end # module