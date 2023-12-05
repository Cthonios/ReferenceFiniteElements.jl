# type used to exploit multiple dispatch
"""
Type to define new element shapes
"""
abstract type ReferenceFEType{N, D, Q} end

"""
Returns the quadrature degree of a ReferenceFEType
"""
degree(e::ReferenceFEType) = e.degree

"""
Returns the number of nodes for a ReferenceFEType
"""
num_nodes(::ReferenceFEType{N, D, Q}) where {N, D, Q} = N

"""
Returns the number of dimensions for a ReferenceFEType
"""
num_dimensions(::ReferenceFEType{N, D, Q}) where {N, D, Q} = D

"""
Returns the number of quadrature porints for a ReferenceFEType
"""
num_q_points(::ReferenceFEType{N, D, Q}) where {N, D, Q} = Q
