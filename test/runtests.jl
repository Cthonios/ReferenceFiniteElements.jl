using Aqua
using Distributions
using LinearAlgebra
using ReferenceFiniteElements
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

function is_inside_hex(point)
  x_cond = (point[1] >= -1.) && (point[1] <= 1.)
  y_cond = (point[2] >= -1.) && (point[2] <= 1.)
  z_cond = (point[3] >= -1.) && (point[3] <= 1.)
  return x_cond && y_cond && z_cond
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
  return points
end

function polyval2d(x::T, y::T, C::Matrix{<:AbstractFloat}) where T
  val = 0.0
  for i in 1:size(C, 1)
    for j in 1:size(C, 2)
      val = val + C[i, j] * x^(i - 1) * y^(j - 1)
    end
  end
  val
end

function dpolyval2d(x::T, y::T, C::Matrix{<:AbstractFloat}, direction::Int) where T
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

function partition_of_unity_shape_function_values_test(el, q_degree, int_type, float_type)
  re = ReferenceFE(el, q_degree, int_type, float_type)
  sums = map(N -> sum(N), re.Ns)
  @test sums ≈ ones(Float64, size(sums))
end

function partition_of_unity_shape_function_gradients_test(el, q_degree, int_type, float_type)
  re = ReferenceFE(el, q_degree, int_type, float_type)
  sums = map(∇N -> sum(∇N), re.∇N_ξs)
  for ∇N in sums
    for i in size(∇N, 1)
      if float_type == Float32
        # @show re.∇N_ξs
        # @show ∇N[i]
        @test ∇N[i] < 2.5e-7
      elseif float_type == Float64
        # @show re.∇N_ξs
        # @show ∇N[i]
        @test ∇N[i] < 1e-13
      end
    end
  end
end

function kronecker_delta_property(el, q_degree, Itype, Rtype)
  e = ReferenceFE(el, q_degree, Itype, Rtype)
  Ns = ReferenceFiniteElements.shape_function_values.((el,), eachcol(e.nodal_coordinates))
  Ns = mapreduce(permutedims, vcat, Ns)
  @test Ns ≈ I
end

function partition_of_unity_tests(el, q_degree, int_type, float_type)
  partition_of_unity_shape_function_values_test(el, q_degree, int_type, float_type)
  partition_of_unity_shape_function_gradients_test(el, q_degree, int_type, float_type)
end

function common_test_sets(el, q_degrees, int_types, float_types)
  for q_degree in q_degrees
    for int_type in int_types
      for float_type in float_types
        @testset ExtendedTestSet "$(typeof(el)), q_degree = $q_degree - Quadrature weight positivity test" begin
          if typeof(el) == Tet4
            continue
          end
          e = ReferenceFE(el, q_degree, int_type, float_type)
          for w in e.ws
            @test w > 0.
          end
        end

        @testset ExtendedTestSet "$(typeof(el)), q_degree = $q_degree - sum of quadrature points test" begin
          e = ReferenceFE(el, q_degree, int_type, float_type)
          if typeof(el) <: ReferenceFiniteElements.HexUnion
            @test e.ws |> sum ≈ 8.0
          elseif typeof(el) <: ReferenceFiniteElements.QuadUnion
            @test e.ws |> sum ≈ 4.0
          elseif typeof(el) <: ReferenceFiniteElements.TetUnion
            @test e.ws |> sum ≈ 1. / 6.
          elseif typeof(el) <: ReferenceFiniteElements.TriUnion
            @test e.ws |> sum ≈ 1. / 2.
          else
            @test false
          end
        end

        @testset ExtendedTestSet "$(typeof(el)), q_degree = $q_degree - triangle exactness" begin
          if typeof(el) <: ReferenceFiniteElements.TriUnion
            # if typeof(el) <: Tri3
            #   continue
            # end
          else
            continue
          end
          
          for degree in [1]
            e = ReferenceFE(el, q_degree)
            for i in 1:degree
              for j in 1:degree - i
                quad_answer = map((ξ, w) -> w * ξ[1]^i * ξ[2]^j, eachcol(e.ξs), e.ws) |> sum
                exact = integrate_2d_monomial_on_triangle(i, j)
                @test quad_answer ≈ exact
              end
            end
          end
        end

        @testset ExtendedTestSet "$(typeof(el)), q_degree = $q_degree - quadrature points inside element" begin
          e = ReferenceFE(el, q_degree, int_type, float_type)
          if typeof(el) <: ReferenceFiniteElements.HexUnion
            test_func = is_inside_hex
          elseif typeof(el) <: ReferenceFiniteElements.QuadUnion
            test_func = is_inside_quad
          elseif typeof(el) <: ReferenceFiniteElements.TetUnion
            test_func = is_inside_tet
          elseif typeof(el) <: ReferenceFiniteElements.TriUnion
            test_func = is_inside_triangle
          else
            @test false
          end

          for ξ in eachcol(e.ξs)
            @test test_func(ξ)
          end
        end

        @testset ExtendedTestSet "$(typeof(el)), q_degree = $q_degree - partition of unity tests" begin
          partition_of_unity_shape_function_values_test(el, q_degree, int_type, float_type)
          partition_of_unity_shape_function_gradients_test(el, q_degree, int_type, float_type)
        end

        @testset ExtendedTestSet "$(typeof(el)), q_degree = $q_degree - kronecker delta property tests" begin
          kronecker_delta_property(el, q_degree, int_type, float_type)
        end
      end
    end
  end
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
