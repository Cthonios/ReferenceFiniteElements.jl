using Aqua
using Distributions
using LinearAlgebra
using ReferenceFiniteElements
using StructArrays
using Test
using TestSetExtensions

function is_inside_unit_interval(point)
  # test is different in optimism from [0., 1.] not [-1., 1.]
  return (point >= 0.) && (point <= 1.)
end

function integrate_1d_monomial_on_line(n)
  return 1. / (n + 1)
end

function map_affine_1d(ξ, end_points)
  return (1. - ξ) * end_points[1] + ξ * end_points[2]
end

function map_1d_jac(end_points)
  return end_points[2] - end_points[1]
end

function is_inside_quad(point)
  x_cond = (point[1] >= -1.) && (point[1] <= 1.)
  y_cond = (point[2] >= -1.) && (point[2] <= 1.)
  return x_cond && y_cond
end

function is_inside_triangle(point)
  x_cond = (point[1] >= 0.) && (point[1] <= 1.)
  y_cond = (point[2] >= 0.) && (point[2] <= 1. - point[1])
  return x_cond && y_cond
end

function is_inside_tet(point)
  x_cond = (point[1] >= 0.) && (point[1] <= 1.)
  y_cond = (point[2] >= 0.) && (point[2] <= 1. - point[1])
  z_cond = (point[3] >= 0.) && (point[3] <= 1. - point[2])
  return x_cond && y_cond && z_cond
end

function integrate_2d_monomial_on_triangle(n, m)
  p = n + m
  return 1. / ((p + 2) * (p + 1) * binomial(p, n))
end

function generate_random_points_in_triangle(n_points::Int)
  x = rand(Uniform(0.0, 1.0), n_points)
  y = zeros(Float64, n_points)
  for i in 1:n_points
    y[i] = rand(Uniform(0.0, 1.0 - x[i]))
  end
  points = hcat(x, y)' |> collect
  # @show points
  return points
end

function polyval2d(x::T, y::T, C::Matrix{T}) where T
  val = 0.0
  for i in 1:size(C, 1)
    for j in 1:size(C, 2)
      val = val + C[i, j] * x^(i - 1) * y^(j - 1)
    end
  end
  val
end

function dpolyval2d(x::T, y::T, C::Matrix{T}, direction::Int) where T
  val = 0.0
  if direction == 1
    for i in 2:size(C, 1)
      for j in 1:size(C, 2)
        val = val + (i - 1) * C[i, j] * x^(i - 2) * y^(j - 1)
      end
    end
  elseif direction == 2
    for i in 1:size(C, 1)
      for j in 2:size(C, 2)
        val = val + (j - 1) * C[i, j] * x^(i - 1) * y^(j - 2)
      end
    end
  end
  val
end

function partition_of_unity_shape_function_values_int_test(el, q_degree)
  shape_functions = ShapeFunctions(el, q_degree)
  sums = map(N -> sum(N), shape_functions.Ns)
  @test sums ≈ ones(Float64, size(sums))
end

function partition_of_unity_shape_function_gradients_int_test(el, q_degree)
  shape_functions = ShapeFunctions(el, q_degree)
  sums = map(∇N -> sum(∇N, dims=1), shape_functions.∇N_ξs)
  for ∇N in sums
    for i in size(∇N, 1)
      # @show ∇N[i]
      @test ∇N[i] < 1e-13
    end
  end
end

function kronecker_delta_property(el, q_degree)
  e = ReferenceFEStencil(el, q_degree)
  Ns = ReferenceFiniteElements.shape_function_values_int.((el,), eachcol(e.coordinates))
  Ns = mapreduce(permutedims, vcat, Ns)
  @test Ns ≈ I
end

function partition_of_unity_tests(el, q_degree)
  partition_of_unity_shape_function_values_int_test(el, q_degree)
  partition_of_unity_shape_function_gradients_int_test(el, q_degree)
end

@includetests ARGS

# Aqua.test_all(ReferenceFiniteElements) # getting weird type ambiguity from StructArrays
Aqua.test_ambiguities(ReferenceFiniteElements)
Aqua.test_unbound_args(ReferenceFiniteElements)
Aqua.test_undefined_exports(ReferenceFiniteElements)
Aqua.test_piracy(ReferenceFiniteElements)
Aqua.test_project_extras(ReferenceFiniteElements)
Aqua.test_stale_deps(ReferenceFiniteElements)
Aqua.test_deps_compat(ReferenceFiniteElements)
Aqua.test_project_toml_formatting(ReferenceFiniteElements)
