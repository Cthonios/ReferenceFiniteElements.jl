using Adapt
using Aqua
using CUDA
using Distributions
using Exodus
using LinearAlgebra
using ReferenceFiniteElements
using StaticArrays
using Test
using TestSetExtensions

function is_inside_element(::ReferenceFiniteElements.AbstractEdge, point)
  return (point[1] >= -1.) && (point[1] <= 1.)
end

function is_inside_surface_element(::ReferenceFiniteElements.AbstractEdge, point)
  return (point[1] >= 0.) && (point[1] <= 1.)
end

function is_inside_element(::ReferenceFiniteElements.AbstractHex, point)
  return (point[1] >= -1.) && (point[1] <= 1.) &&
         (point[2] >= -1.) && (point[2] <= 1.) &&
         (point[3] >= -1.) && (point[3] <= 1.)
end

function is_inside_element(::ReferenceFiniteElements.AbstractQuad, point)
  return (point[1] >= -1.) && (point[1] <= 1.) &&
         (point[2] >= -1.) && (point[2] <= 1.)
end

function is_inside_element(::ReferenceFiniteElements.AbstractTri, point)
  return (point[1] >= 0.) && (point[1] <= 1.) &&
         (point[2] >= 0.) && (point[2] <= 1. - point[1])
end

function is_inside_element(::ReferenceFiniteElements.AbstractTet, point)
  return (point[1] >= 0.) && (point[1] <= 1.) &&
         (point[2] >= 0.) && (point[2] <= 1. - point[1]) && 
         (point[3] >= 0.) && (point[3] <= 1. - point[2])
end


q_weight_sum(::ReferenceFiniteElements.AbstractEdge) = 2.
q_weight_sum(::ReferenceFiniteElements.AbstractHex) = 8.
q_weight_sum(::ReferenceFiniteElements.AbstractQuad) = 4.
q_weight_sum(::ReferenceFiniteElements.AbstractTet) = 1. / 6.
q_weight_sum(::ReferenceFiniteElements.AbstractTri) = 0.5
q_weight_sum(::ReferenceFiniteElements.AbstractVertex) = 1.

function test_q_points_inside_element(re) 
  @test all(is_inside_element.((re.element,), quadrature_points(re)))
  for n in 1:num_quadrature_points(re)
    @test is_inside_element(re.element, quadrature_point(re, n))
  end
  @test all(is_inside_element.((re.element,), surface_quadrature_points(re)))
  for f in 1:num_faces(re.element)
    for n in 1:num_quadrature_points(re.surface_element)
      if typeof(re.element) <: ReferenceFiniteElements.AbstractTri
        @test is_inside_element(ReferenceFiniteElements.surface_element(re.element), 2. * surface_quadrature_point(re, n, f) .- 1.)
      else
        @test is_inside_element(ReferenceFiniteElements.surface_element(re.element), surface_quadrature_point(re, n, f))
      end
    end
  end
end

function test_q_weight_positivity(re) 
  ws = quadrature_weights(re)
  @test all(ws .> zero(eltype(ws)))
  for n in 1:num_quadrature_points(re)
    @test quadrature_weight(re, n) > zero(eltype(ws))
  end
  ws = surface_quadrature_weights(re)
  @test all(ws .> zero(eltype(ws)))
  for f in 1:num_faces(re.element)
    for n in 1:num_quadrature_points(re.surface_element)
      @test surface_quadrature_weight(re, n, f) > zero(eltype(ws))
    end
  end
end

function test_q_weight_sum(re) 
  @test q_weight_sum(re.element) ≈ sum(quadrature_weights(re))
  for f in 1:num_faces(re.element)
    if typeof(re.element) <: ReferenceFiniteElements.AbstractTri
      @test sum(surface_quadrature_weights(re)[:, f]) ≈ q_weight_sum(ReferenceFiniteElements.surface_element(re.element)) / 2.
    else
      @test sum(surface_quadrature_weights(re)[:, f]) ≈ q_weight_sum(ReferenceFiniteElements.surface_element(re.element))
    end
  end
end 

# TODO need a test on exact integration for quadrature rules
function test_partition_of_unity_on_values(re) 
  Ns = shape_function_values(re)
  @test all(sum.(Ns) .≈ one(eltype(Ns[1])))
  for n in 1:num_quadrature_points(re)
    N = shape_function_value(re, n)
    @test sum(N) ≈ one(eltype(N))
  end

  Ns = surface_shape_function_values(re)
  @test all(sum.(Ns) .≈ one(eltype(Ns[1, 1])))
  for f in 1:num_faces(re.element)
    for n in 1:num_quadrature_points(ReferenceFiniteElements.surface_element(re.element))
      N = surface_shape_function_value(re, n, f)
      @test sum(N) ≈ one(eltype(N))
    end
  end
end

function test_partition_of_unity_on_gradients(re) 
  for ∇N in shape_function_gradients(re)
    @test isapprox(sum(∇N), 0, atol=5e-12)
  end
  for n in 1:num_quadrature_points(re)
    ∇N = shape_function_gradient(re, n)
    @test isapprox(sum(∇N), 0, atol=5e-12)
  end
  for ∇N in surface_shape_function_gradients(re)
    @test isapprox(sum(∇N), 0, atol=5e-12)
  end
  for f in 1:num_faces(re.element)
    for n in 1:num_quadrature_points(ReferenceFiniteElements.surface_element(re.element))
      ∇N = surface_shape_function_gradient(re, n, f)
      @test isapprox(sum(∇N), 0, atol=5e-12)
    end
  end
end

function test_partition_of_unity_on_hessians(re)
  for ∇∇N in shape_function_hessians(re)
    @test isapprox(sum(∇∇N), 0, atol=5e-11)
  end
  for n in 1:num_quadrature_points(re)
    ∇∇N = shape_function_hessian(re, n)
    @test isapprox(sum(∇∇N), 0, atol=5e-11)
  end
  for ∇∇N in surface_shape_function_hessians(re)
    @test isapprox(sum(∇∇N), 0, atol=5e-11)
  end
  for f in 1:num_faces(re.element)
    for n in 1:num_quadrature_points(ReferenceFiniteElements.surface_element(re.element))
      ∇∇N = surface_shape_function_hessian(re, n, f)
      @test isapprox(sum(∇∇N), 0, atol=5e-11)
    end
  end
end

function test_kronecker_delta_property(re) 
  if ReferenceFiniteElements.polynomial_degree(re.element) == 0
    @test isapprox(hcat(map(x -> 
      shape_function_value(re.element, re.Xs, x, re.backend), re.Xs)...), 
      fill(1., (ReferenceFiniteElements.num_shape_functions(re.element), length(re.Xs))), 
      atol=5e-14
    )
  else
    @test isapprox(hcat(map(x -> 
      shape_function_value(re.element, re.Xs, x, re.backend), re.Xs)...), 
      I, 
      atol=5e-14
    )
  end
  # TODO add test on surface shape functions
  # will need to first calculate shape function then 
  # index based on the edges, since non-edge dofs will be zero 
end

# extension tests
include("TestAdaptExt.jl")
include("TestExodusExt.jl")

# main package tests
# el_types = [Edge2, Edge2, Edge3, Quad4, Quad4]
# q_orders = [1, 2, 2, 1, 2]
el_types = [
  (Edge0, 1),
  (Edge2, 1),
  (Edge2, 2),
  (Edge3, 1),
  (Edge3, 2),
  (Hex0, 1),
  (Hex8, 1),
  (Hex8, 2),
  (Quad0, 1),
  (Quad4, 1),
  (Quad4, 2),
  (Quad9, 1),
  (Quad9, 2),
  (Tet0, 1),
  (Tet4, 1),
  (Tet4, 2),
  (Tet10, 1),
  (Tet10, 2),
  # (Tri0, 1),
  (Tri3, 1),
  (Tri3, 2),
  (Tri6, 1),
  (Tri6, 2)
]
# for (el_type, q_order) in zip(el_types, q_orders)
for el_type in el_types
  el_type, q_order = el_type
  @testset "$el_type - q order = $q_order Tests" begin
    re = ReferenceFE{Int64, Float64, SArray}(el_type{Lagrange, q_order}())
    test_q_points_inside_element(re)
    if el_type <: ReferenceFiniteElements.AbstractTet || q_order > 1
      # do nothing here
    else
      test_q_weight_positivity(re)
    end
    test_q_weight_sum(re)
    test_partition_of_unity_on_values(re)
    test_partition_of_unity_on_gradients(re)
    test_partition_of_unity_on_hessians(re)
    test_kronecker_delta_property(re)
  end
end

el_types = [Edge, Quad]
for el_type in el_types
  @testset "$el_type Tests" begin
    qs = 1:10
    for q in qs
      re = ReferenceFE{SArray, el_type, Lagrange, 1, q}()
      test_q_points_inside_element(re)
      test_q_weight_positivity(re)
      test_q_weight_sum(re)
    end

    # just test initializing this
    # re = ReferenceFE{SArray, el_type, Lagrange, 0, 1}()
    # test_partition_of_unity_on_values(re)
    # test_partition_of_unity_on_gradients(re)
    # test_partition_of_unity_on_hessians(re)
    # test_kronecker_delta_property(re)
    
    ps = 1:10
    for p in ps
      if p > 5
        arr_type = Array
      else
        arr_type = SArray
      end
      re = ReferenceFE{arr_type, el_type, Lagrange, p, p}()
      test_partition_of_unity_on_values(re)
      test_partition_of_unity_on_gradients(re)
      test_partition_of_unity_on_hessians(re)
      test_kronecker_delta_property(re)
    end
  end
end

@testset ExtendedTestSet "Aqua Tests" begin
  Aqua.test_all(ReferenceFiniteElements; ambiguities=false)
end

# JET testing
# @testset ExtendedTestSet "JET Tests" begin
#   # invalidations from FastGaussQuadrature so only targeting defined
#   # modules
#   test_package("ReferenceFiniteElements"; target_defined_modules=true)
# end
