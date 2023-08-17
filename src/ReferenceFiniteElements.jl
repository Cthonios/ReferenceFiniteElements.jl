module ReferenceFiniteElements

# element types
# export Edge
export Hex8
export Quad4, Quad9
export Tet4, Tet10
export Tri3, Tri6

# types
export ReferenceFE

# methods
export quadrature_point
export quadrature_points
export quadrature_weight
export quadrature_weights
export shape_function_gradients
export shape_function_hessians
export shape_function_values
export vertex_nodes

# dependencies
using FastGaussQuadrature
using InteractiveUtils
using LinearAlgebra
using Polynomials
using PrecompileTools
using SpecialPolynomials
using StaticArrays
using StructArrays

# type used to exploit multiple dispatch
abstract type ReferenceFEType{N, D} end
degree(e::ReferenceFEType) = e.degree
num_nodes(::ReferenceFEType{N, D}) where {N, D} = N
num_dimensions(::ReferenceFEType{N, D}) where {N, D} = D

# type to aid in making ReferenceFE with StructArrays.jl
struct Interpolants{N, D, Ftype, L1, L2}
  ξ::SVector{D, Ftype}
  w::Ftype
  N::SVector{N, Ftype}
  ∇N_ξ::SMatrix{N, D, Ftype, L1}
  ∇∇N_ξ::SArray{Tuple{N, D, D}, Ftype, 3, L2}
end

function Interpolants(
  e::ReferenceFEType{N, D}, ::Type{Ftype} = Float64
) where {N, D, Ftype}

  ξs_temp, ws = quadrature_points_and_weights(e, Ftype)
  ξs = reinterpret(SVector{D, Ftype}, vec(ξs_temp))
  Ns = Vector{SVector{N, Ftype}}(undef, length(ξs))
  ∇N_ξs = Vector{SMatrix{N, D, Ftype, N * D}}(undef, length(ξs))
  ∇∇N_ξs = Vector{SArray{Tuple{N, D, D}, Ftype, 3, N * D * D}}(undef, length(ξs))
  for (q, ξ) in enumerate(ξs)
    Ns[q]     = shape_function_values(e, ξ)
    ∇N_ξs[q]  = shape_function_gradients(e, ξ)
    ∇∇N_ξs[q] = shape_function_hessians(e, ξ)
  end
  s = StructArray{Interpolants{N, D, Ftype, N * D, N * D * D}}((
    ξ=ξs, w=ws,
    N=Ns, ∇N_ξ=∇N_ξs, ∇∇N_ξ=∇∇N_ξs
  ))
  return s
end

# main type of the package
struct ReferenceFE{Itype, N, D, Ftype, L1, L2}
  nodal_coordinates::VecOrMat{Ftype}
  face_nodes::Matrix{Itype}
  interior_nodes::Vector{Itype}
  # figure out a way to fix below, that's dumb
  interpolants::StructVector{
    Interpolants{N, D, Ftype, L1, L2}, 
    NamedTuple{
      (:ξ, :w, :N, :∇N_ξ, :∇∇N_ξ), 
      Tuple{
        Vector{SVector{D, Ftype}}, 
        Vector{Ftype}, 
        Vector{SVector{N, Ftype}}, 
        Vector{SMatrix{N, D, Ftype, L1}},
        Vector{SArray{Tuple{N, D, D}, Ftype, 3, L2}}
      }
    }, 
    Int64
  }
end

function ReferenceFE(
  e::ReferenceFEType{N, D},
  ::Type{Itype} = Int64, ::Type{Ftype} = Float64
) where {N, D, Itype, Ftype}

  nodal_coordinates, face_nodes, interior_nodes = element_stencil(e, Itype, Ftype)
  interps = Interpolants(e, Ftype)

  return ReferenceFE{Itype, N, D, Ftype, N * D, N * D *D}(
    nodal_coordinates, face_nodes, interior_nodes,
    interps
  )
end

function Base.show(io::IO, e::ReferenceFE)
  print(io, "Element type             = $(typeof(e))\n\n")
  print(io, "Nodal coordinates        = \n")
  display(e.nodal_coordinates)
  print(io, "\n")
  print(io, "Face nodes               = \n")
  display(e.face_nodes)
  print(io, "\n")
  print(io, "Interior nodes           = \n")
  display(e.interior_nodes)
  print(io, "\n")
  print(io, "Shape function values    = \n")
  for n in axes(e.interpolants, 1)
    display(shape_function_values(e, n))
  end
  print(io, "\n")
  print(io, "Shape function gradients = \n")
  for n in axes(e.interpolants, 1)
    display(shape_function_gradients(e, n))
  end
  print(io, "\n")
  print(io, "Shape function hessians  = \n")
  for n in axes(e.interpolants, 1)
    display(shape_function_hessians(e, n))
  end
  print(io, "\n")
end

quadrature_point(e::ReferenceFE, q::Integer) = LazyRow(getfield(e, :interpolants), q).ξ
quadrature_points(e::ReferenceFE) = getfield(e, :interpolants).ξ
quadrature_weight(e::ReferenceFE, q::Integer) = LazyRow(getfield(e, :interpolants), q).w
quadrature_weights(e::ReferenceFE) = getfield(e, :interpolants).w
shape_function_values(e::ReferenceFE) = getfield(e, :interpolants).N
shape_function_values(e::ReferenceFE, i::Integer) = LazyRow(getfield(e, :interpolants), i).N
shape_function_gradients(e::ReferenceFE) = getfield(e, :interpolants).∇N_ξ
shape_function_gradients(e::ReferenceFE, i::Integer) = LazyRow(getfield(e, :interpolants), i).∇N_ξ
shape_function_hessians(e::ReferenceFE) = getfield(e, :interpolants).∇∇N_ξ
shape_function_hessians(e::ReferenceFE, i::Integer) = LazyRow(getfield(e, :interpolants), i).∇∇N_ξ
vertex_nodes(::ReferenceFE{Itype, N, D, Ftype, L1, L2}) where {Itype, N, D, Ftype, L1, L2} = 1:N

# implementations of things common across multiple element types
include("implementations/HexCommon.jl")
include("implementations/QuadCommon.jl")
include("implementations/TetCommon.jl")
include("implementations/TriCommon.jl")

# implementations of things specific to element types
# include("implementations/Edge.jl")
include("implementations/Hex8.jl")
include("implementations/Quad4.jl")
include("implementations/Quad9.jl")
include("implementations/Tet4.jl")
include("implementations/Tet10.jl")
include("implementations/Tri3.jl")
include("implementations/Tri6.jl")

# include("implementations/SimplexTri.jl")

# precompilation
@setup_workload begin
  @compile_workload begin
    for degree in [1, 2, 3, 4, 5, 6]
      ReferenceFE(Hex8(degree))
      ReferenceFE(Quad4(degree))
      ReferenceFE(Quad9(degree))
      # ReferenceFE(Tet4(degree))
      # ReferenceFE(Tet10(degree))
      ReferenceFE(Tri3(degree))
      ReferenceFE(Tri6(degree))
    end

    for degree in [1, 2]
      ReferenceFE(Tet4(degree))
    end
  end
end

end # module ReferenceFEs
