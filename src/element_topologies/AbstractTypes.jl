"""
$(TYPEDEF)
"""
abstract type AbstractElementTopology{D, V, E, F} end
"""
$(TYPEDSIGNATURES)
"""
ndimension(::AbstractElementTopology{D, V, E, F}) where {D, V, E, F} = D
"""
$(TYPEDSIGNATURES)
"""
nvertices(::AbstractElementTopology{D, V, E, F}) where {D, V, E, F} = V
"""
$(TYPEDSIGNATURES)
"""
nedges(::AbstractElementTopology{D, V, E, F}) where {D, V, E, F} = E
"""
$(TYPEDSIGNATURES)
"""
nfaces(::AbstractElementTopology{D, V, E, F}) where {D, V, E, F} = F

# other methods to define for subtypes
function ninteriorvertices end
function nverticesperedge end
function nverticesperface end
# function surfacelement end

# function surface_element(::Abstract{D, V, E, F})

"""
$(TYPEDEF)
"""
abstract type AbstractVertex <: AbstractElementTopology{0, 1, 0, 0} end
surface_element(::AbstractVertex) = nothing

"""
$(TYPEDEF)
"""
abstract type AbstractEdge{V} <: AbstractElementTopology{1, V, 1, 0} end
faces(::AbstractEdge) = Int[Int[]]

"""
$(TYPEDEF)
"""
abstract type AbstractFace{V, E} <: AbstractElementTopology{2, V, E, 1} end

"""
$(TYPEDEF)
"""
abstract type AbstractQuad{V} <: AbstractElementTopology{2, V, 4, 1} end

"""
$(TYPEDEF)
"""
abstract type AbstractTri{V} <: AbstractElementTopology{2, V, 3, 1} end

"""
$(TYPEDEF)
"""
abstract type AbstractVolume{V, E, F} <: AbstractElementTopology{3, V, E, F} end

"""
$(TYPEDEF)
"""
abstract type AbstractHex{V} <: AbstractElementTopology{3, V, 12, 6} end

"""
$(TYPEDEF)
"""
abstract type AbstractTet{V} <: AbstractElementTopology{3, V, 6, 4} end
