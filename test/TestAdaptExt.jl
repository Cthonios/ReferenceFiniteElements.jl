function test_adapt_cuda(re_dev)
  @test device(re_dev.Xs) |> typeof <: CuDevice
  @test device(re_dev.edge_nodes) |> typeof <: CuDevice
  @test device(re_dev.face_nodes) |> typeof <: CuDevice
  @test device(re_dev.interior_nodes) |> typeof <: CuDevice
  @test device(re_dev.cell_interps.vals.N) |> typeof <: CuDevice
  @test device(re_dev.cell_interps.vals.w) |> typeof <: CuDevice
  @test device(re_dev.cell_interps.vals.ξ) |> typeof <: CuDevice
  @test device(re_dev.cell_interps.vals.∇N_ξ) |> typeof <: CuDevice
  @test device(re_dev.cell_interps.vals.∇∇N_ξ) |> typeof <: CuDevice
end

function test_adapt_rocm(re_dev)
  # @test device(re_dev.Xs) |> typeof <: ROCDeviceArray
  # @test device(re_dev.edge_nodes) |> typeof <: ROCDeviceArray
  # @test device(re_dev.face_nodes) |> typeof <: ROCDeviceArray
  # @test device(re_dev.interior_nodes) |> typeof <: ROCDeviceArray
  # @test device(re_dev.cell_interps.vals.N) |> typeof <: ROCDeviceArray
  # @test device(re_dev.cell_interps.vals.w) |> typeof <: ROCDeviceArray
  # @test device(re_dev.cell_interps.vals.ξ) |> typeof <: ROCDeviceArray
  # @test device(re_dev.cell_interps.vals.∇N_ξ) |> typeof <: ROCDeviceArray
  # @test device(re_dev.cell_interps.vals.∇∇N_ξ) |> typeof <: ROCDeviceArray
end

@testset ExtendedTestSet "AdaptExt - CPU" begin
  re = ReferenceFE(Quad4{Lagrange, 2}())
  re = Adapt.adapt_structure(Array, re)
end

@testset ExtendedTestSet "AdaptExt - AMDGPU" begin
  if AMDGPU.has_rocm_gpu()
    re = ReferenceFE(Quad4{Lagrange, 2}())
    re = Adapt.adapt_structure(ROCArray, re)
    display(re)
    test_adapt_rocm(re)
  end
end

@testset ExtendedTestSet "AdaptExt - CUDA" begin
  if CUDA.has_cuda()
    re = ReferenceFE(Quad4{Lagrange, 2}())
    re = Adapt.adapt_structure(CuArray, re)
    test_adapt_cuda(re)
  end
end
