module ReferenceFiniteElements

# element types
# export Edge
export Hex8,
       Quad4, Quad9,
       Tet4, Tet10,
       Tri3, Tri6

# types
export ReferenceFE

# methods.
export num_dimensions,
       num_nodes,
       quadrature_point,
       quadrature_points,
       quadrature_weight,
       quadrature_weights,
       shape_function_gradients,
       shape_function_hessians,
       shape_function_values,
       vertex_nodes

# dependencies
using DocStringExtensions,
      FastGaussQuadrature,
      LinearAlgebra,
      Polynomials,
      PrecompileTools,
      SpecialPolynomials,
      StaticArrays,
      StructArrays

# for docs
@template (FUNCTIONS, METHODS, MACROS) =
"""
$(TYPEDSIGNATURES)
$(DOCSTRING)
$(METHODLIST)
"""

@template (TYPES) =
"""
$(TYPEDFIELDS)
$(DOCSTRING)
"""

# Types
include("ReferenceFETypes.jl")
include("Interpolants.jl")
include("ReferenceFEs.jl")
include("Utils.jl")

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
    for int_type in [Int32, Int64]
      for float_type in [Float32, Float64]
        for array_type in [SArray, MArray]
          for degree in [1, 2, 3, 4, 5, 6]
            ReferenceFE(Hex8(degree); int_type=int_type, float_type=float_type, array_type=array_type)
            ReferenceFE(Quad4(degree); int_type=int_type, float_type=float_type, array_type=array_type)
            ReferenceFE(Quad9(degree); int_type=int_type, float_type=float_type, array_type=array_type)
            ReferenceFE(Tri3(degree); int_type=int_type, float_type=float_type, array_type=array_type)
            ReferenceFE(Tri6(degree); int_type=int_type, float_type=float_type, array_type=array_type)
          end

          for degree in [1, 2]
            ReferenceFE(Tet4(degree); int_type=int_type, float_type=float_type, array_type=array_type)
            ReferenceFE(Tet10(degree); int_type=int_type, float_type=float_type, array_type=array_type)
          end
        end
      end
    end
  end
end

end # module ReferenceFEs
