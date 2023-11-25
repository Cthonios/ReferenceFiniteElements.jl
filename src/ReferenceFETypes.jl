# type used to exploit multiple dispatch
"""
Type to define new element shapes
"""
abstract type ReferenceFEType{N, D} end

"""
Returns the quadrature degree of a ReferenceFEType
"""
degree(e::ReferenceFEType) = e.degree

"""
Returns the number of nodes for a ReferenceFEType
"""
num_nodes(::ReferenceFEType{N, D}) where {N, D} = N

"""
Returns the number of dimensions for a ReferenceFEType
"""
num_dimensions(::ReferenceFEType{N, D}) where {N, D} = D