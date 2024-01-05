using Adapt
using Aqua
using CUDA
using Distributions
using Exodus
using ForwardDiff
using JET
using KernelAbstractions
using LinearAlgebra
using ReferenceFiniteElements
using StaticArrays
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

@kernel function sum_shape_function_values_kernel!(sums, shapes)
  I = @index(Global)
  sums[I] = sum(shapes[I])
end

function partition_of_unity_shape_function_values_test(re, backend)
  shapes = shape_function_values(re)
  sums = KernelAbstractions.zeros(backend, ReferenceFiniteElements.float_type(re), length(shapes))
  kernel = sum_shape_function_values_kernel!(backend)
  kernel(sums, shapes, ndrange=length(shapes))
  sums_cpu = convert(Vector, sums)
  @test sums_cpu ≈ ones(eltype(sums_cpu), length(sums_cpu))
end

@kernel function sum_shape_function_gradients_kernel!(sums, shapes)
  I = @index(Global)
  sums[I] = sum(shapes[I])
end

function partition_of_unity_shape_function_gradients_test(re, backend)
  shapes = shape_function_gradients(re)
  sums = KernelAbstractions.zeros(backend, ReferenceFiniteElements.float_type(re), length(shapes))
  kernel = sum_shape_function_gradients_kernel!(backend)
  kernel(sums, shapes, ndrange=length(sums))
  sums_cpu = convert(Vector, sums)
  for ∇N in sums_cpu
    for i in size(∇N, 1)
      if eltype(sums_cpu) == Float32
        @test ∇N[i] < 2.5e-7
      elseif eltype(sums_cpu) == Float64
        @test ∇N[i] < 1e-13
      end
    end
  end
end

@kernel function kronecker_delta_property_kernel!(Ns, re, coords, mvec=false)
  I = @index(Global)
  # note only SVector makes sense here. MVector is no bueno in kernels
  if mvec
    type = MVector
  else
    type = SVector
  end
  Ns[I] = ReferenceFiniteElements.shape_function_values(typeof(re.ref_fe_type)(num_q_points(re)), type, coords[I])
end

function kronecker_delta_property(re, backend, arr_type)
  n_dim = size(re.nodal_coordinates, 1)
  coords = reinterpret(SVector{n_dim, ReferenceFiniteElements.float_type(re)}, re.nodal_coordinates)
  Ns = KernelAbstractions.zeros(backend, arr_type{num_nodes_per_element(re), eltype(re.nodal_coordinates)}, length(coords))
  kernel = kronecker_delta_property_kernel!(backend)
  kernel(Ns, re, coords, ndrange=length(Ns))
  Ns = convert(Vector, Ns)
  Ns = mapreduce(permutedims, vcat, Ns)
  @test Ns ≈ I
end

# just so Enzyme can autodiff below. This should probably be in StaticArrays.jl
Base.one(::Type{SVector{N, T}}) where {N, T} = ones(SVector{N, T})

@kernel function test_gradients_kernel!(re, ξs, ∇Ns_fd)
  I = @index(Global)

  # if mvec
  #   type = MVector
  # else
  #   type = SVector
  # end
  type = SVector
  # ∇Ns_fd[I] = Zygote.jacobian(x -> shape_function_values(re.ref_fe_type, type, x), ξs[I])[1]
  # ∇Ns_fd[I] = AD.jacobian(AD.ZygoteBackend(), x -> shape_function_values(re.ref_fe_type, type, x), ξs[I])[1]
  # ∇Ns_fd[I] = Enzyme.jacobian(Reverse, x -> shape_function_values(re.ref_fe_type, type, x), ξs[I], Val{8})
  # temp = autodiff(Reverse, shape_function_values, re.ref_fe_type, type, Active(ξs[I]))
  temp = Enzyme.jacobian(Reverse, x -> shape_function_values(re.ref_fe_type, type, x), ξs[I], Val{3}())
  @show temp
  # ∇Ns_fd[I] = temp[3]  
end

# function test_gradients(re, backend)
function test_gradients(re)
  # TODO hardcoding to SVector for now
  # ξs = quadrature_points(re)
  # ∇Ns_an = shape_function_gradients(re)
  # # ∇Ns_fd = Vector{eltype(∇Ns_an)}(undef, length(∇Ns_an))
  # ∇Ns_fd = KernelAbstractions.zeros(backend, eltype(∇Ns_an), length(∇Ns_an))
  # # ∇Ns_fd = Adapt.adapt_structure()
  # kernel = test_gradients_kernel!(backend)
  # kernel(re, ξs, ∇Ns_fd, ndrange=length(∇Ns_fd))
  # ∇Ns_an = convert(Vector, ∇Ns_an)
  # ∇Ns_fd = convert(Vector, ∇Ns_fd)
  # for q in axes(∇Ns_an, 1)
  #   @test ∇Ns_an[q] ≈ ∇Ns_fd[q]
  # end
  type = SVector
  for (q, ξ) in enumerate(quadrature_points(re))
    ∇Ns_an = shape_function_gradients(re, q)
    ∇Ns_fd = ForwardDiff.jacobian(x -> shape_function_values(re.ref_fe_type, type, x), ξ)
    @test ∇Ns_fd ≈ ∇Ns_an
  end
end

function test_hessians(re)
  # TODO hardcoding to SMatrix for now
  for (q, ξ) in enumerate(quadrature_points(re))
    ∇∇Ns_an = shape_function_hessians(re, q)
    ∇∇Ns_fd = reshape(ForwardDiff.jacobian(x -> shape_function_gradients(re.ref_fe_type, SMatrix, x), ξ), size(∇∇Ns_an))
    @test ∇∇Ns_an ≈ ∇∇Ns_fd
  end
end

function test_show(re)
  @show re
end

function test_adapt_cuda(re_cuda)
  @test device(re_cuda.nodal_coordinates) |> typeof <: CuDevice
  @test device(re_cuda.face_nodes) |> typeof <: CuDevice
  @test device(re_cuda.interior_nodes) |> typeof <: CuDevice
  @test device(re_cuda.interpolants.N) |> typeof <: CuDevice
  @test device(re_cuda.interpolants.w) |> typeof <: CuDevice
  @test device(re_cuda.interpolants.ξ) |> typeof <: CuDevice
  @test device(re_cuda.interpolants.∇N_ξ) |> typeof <: CuDevice
  @test device(re_cuda.interpolants.∇∇N_ξ) |> typeof <: CuDevice
end

# @kernel function test_quadrature_weight_positivity_kernel(ws)
#   I = @index(Global)
#   # @test ws[I] > zero(eltype(ws))
# end

function test_quadrature_weight_positivity(re, backend)
  ws = quadrature_weights(re)
  ws = convert(Vector, ws)
  for w in ws
    @test w > zero(ReferenceFiniteElements.float_type(re))
  end
  # for w in quadrature_weights(re)
    # @test w > 0.
  # end
  # kernel = test_quadrature_weight_positivity_kernel(backend)
  # kernel(ws, ndrange=length(ws))
end

function common_test_sets_inner(el, q_degree, int_type, float_type, array_type, storage_type, cuda, re)
  if cuda
    re = Adapt.adapt_structure(CuArray, re)
    backend = CUDABackend()
  else
    backend = CPU()
  end

  @testset ExtendedTestSet "$el, $q_degree, $int_type, $float_type, $array_type, $storage_type - Quadrature weight positivity test" begin
    if typeof(re.ref_fe_type) <: Tet4 || typeof(re.ref_fe_type) <: Tet10
      
    else
      # e = ReferenceFE(el(q_degree); int_type, float_type, array_type)
      # for w in quadrature_weights(re)
      #   @test w > 0.
      # end
      test_quadrature_weight_positivity(re, backend)
    end
  end

  @testset ExtendedTestSet "$el, $q_degree, $int_type, $float_type, $array_type, $storage_type - sum of quadrature points test" begin
    # e = ReferenceFE(el(q_degree); int_type, float_type, array_type)
    if typeof(el(q_degree)) <: ReferenceFiniteElements.AbstractHex
      @test quadrature_weights(re) |> sum ≈ 8.0
    elseif typeof(el(q_degree)) <: ReferenceFiniteElements.AbstractQuad
      @test quadrature_weights(re) |> sum ≈ 4.0
    elseif typeof(el(q_degree)) <: ReferenceFiniteElements.AbstractTet
      @test quadrature_weights(re) |> sum ≈ 1. / 6.
    elseif typeof(el(q_degree)) <: ReferenceFiniteElements.AbstractTri
      @test quadrature_weights(re) |> sum ≈ 1. / 2.
    else
      @test false
    end
  end

  @testset ExtendedTestSet "$el, $q_degree, $int_type, $float_type, $array_type, $storage_type - triangle exactness" begin
    if typeof(el) <: ReferenceFiniteElements.AbstractTri
      for degree in [1]
        e = ReferenceFE(el(q_degree); int_type, float_type, array_type)
        for i in 1:degree
          for j in 1:degree - i
            quad_answer = map((ξ, w) -> w * ξ[1]^i * ξ[2]^j, eachcol(e.ξs), e.ws) |> sum
            exact = integrate_2d_monomial_on_triangle(i, j)
            @test quad_answer ≈ exact
          end
        end
      end
    else

    end
  end

  @testset ExtendedTestSet "$el, $q_degree, $int_type, $float_type, $array_type, $storage_type - quadrature points inside element" begin
    e = ReferenceFE(el(q_degree); int_type, float_type, array_type)
    if typeof(el(q_degree)) <: ReferenceFiniteElements.AbstractHex
      test_func = is_inside_hex
    elseif typeof(el(q_degree)) <: ReferenceFiniteElements.AbstractQuad
      test_func = is_inside_quad
    elseif typeof(el(q_degree)) <: ReferenceFiniteElements.AbstractTet
      test_func = is_inside_tet
    elseif typeof(el(q_degree)) <: ReferenceFiniteElements.AbstractTri
      test_func = is_inside_triangle
    else
      @test false
    end

    for ξ in quadrature_points(e)
      @test test_func(ξ)
    end
  end

  @testset ExtendedTestSet "$el, $q_degree, $int_type, $float_type, $array_type, $storage_type - partition of unity tests" begin
    partition_of_unity_shape_function_values_test(re, backend)
    partition_of_unity_shape_function_gradients_test(re, backend)
  end

  @testset ExtendedTestSet "$el, $q_degree, $int_type, $float_type, $array_type, $storage_type - kronecker delta property tests" begin
    if array_type == SArray
      kronecker_delta_property(re, backend, SVector)
    else
      kronecker_delta_property(re, backend, MVector)
    end
  end

  # TODO find an AD engine that actually works well with GPUarrays and staticarrays
  if !cuda
    @testset ExtendedTestSet "$el, $q_degree, $int_type, $float_type, $array_type, $storage_type - shape function gradients" begin
      test_gradients(re)
    end

    @testset ExtendedTestSet "$el, $q_degree, $int_type, $float_type, $array_type, $storage_type - shape function hessians" begin
      test_hessians(re)
    end
  end

  @testset ExtendedTestSet "$el, $q_degree, $int_type, $float_type, $array_type, $storage_type - Base.show" begin
    if !cuda
      test_show(re)
    end
  end

  if cuda
    if array_type <: SArray
      @testset ExtendedTestSet "$el, $q_degree, $int_type, $float_type, $array_type, $storage_type - AdaptExt, CUDA" begin
        test_adapt_cuda(re)
      end
    else
      @info "Skipping mutable static arrays since these are not supported on CUDA"
    end
  end
end

function common_test_sets(el, q_degrees, int_types, float_types, array_types, storage_types; cuda = false)
  for q_degree in q_degrees
    for int_type in int_types
      for float_type in float_types
        for array_type in array_types
          for storage_type in storage_types
            re = ReferenceFE(el(q_degree); int_type, float_type, array_type, storage_type)
            common_test_sets_inner(el, q_degree, int_type, float_type, array_type, storage_type, cuda, re)
          end
        end
      end
    end
  end # q_degree
end # common_test_sets

@includetests ARGS

@testset ExtendedTestSet "Aqua Tests" begin
  # ambiguities due to StructArrays
  Aqua.test_all(ReferenceFiniteElements; ambiguities=false)
end

# JET testing
@testset ExtendedTestSet "JET Tests" begin
  # lots of erros if we don't target only ReferenceFiniteElements
  test_package("ReferenceFiniteElements"; target_defined_modules=true)
end
