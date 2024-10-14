using Aqua
using Distributions
using JET
using LinearAlgebra
using ReferenceFiniteElements
using StaticArrays
using Test
using TestSetExtensions

function is_inside_element(::Edge, point)
  return (point[1] >= -1.) && (point[1] <= 1.)
end

function is_inside_element(::Quad, point)
  return (point[1] >= -1.) && (point[1] <= 1.) &&
         (point[2] >= -1.) && (point[2] <= 1.)
end

q_weight_sum(::Edge) = 2.
q_weight_sum(::Quad) = 4.
q_weight_sum(::Tri) = 0.5
q_weight_sum(::Vertex) = 1.

function test_q_points_inside_element(re) 
  @test all(is_inside_element.((re.element,), quadrature_points(re)))
  for n in 1:num_quadrature_points(re)
    @test is_inside_element(re.element, quadrature_point(re, n))
  end
  @test all(is_inside_element.((re.element,), surface_quadrature_points(re)))
  for f in 1:num_faces(re.element)
    for n in 1:num_quadrature_points(re.surface_element)
      @test is_inside_element(ReferenceFiniteElements.surface_element(re.element), surface_quadrature_point(re, n, f))
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
    @test sum(surface_quadrature_weights(re)[:, f]) ≈ q_weight_sum(ReferenceFiniteElements.surface_element(re.element))
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
    re = ReferenceFE{SArray, el_type, Lagrange, 0, 1}()
    test_partition_of_unity_on_values(re)
    test_partition_of_unity_on_gradients(re)
    test_partition_of_unity_on_hessians(re)
    test_kronecker_delta_property(re)
    
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
@testset ExtendedTestSet "JET Tests" begin
  # invalidations from FastGaussQuadrature so only targeting defined
  # modules
  test_package("ReferenceFiniteElements"; target_defined_modules=true)
end
