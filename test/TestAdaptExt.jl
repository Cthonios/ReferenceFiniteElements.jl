function test_adapt_cuda(re_cuda)
  @test device(re_cuda.Xs) |> typeof <: CuDevice
  @test device(re_cuda.edge_nodes) |> typeof <: CuDevice
  @test device(re_cuda.face_nodes) |> typeof <: CuDevice
  @test device(re_cuda.interior_nodes) |> typeof <: CuDevice
  @test device(re_cuda.cell_interps.vals.N) |> typeof <: CuDevice
  @test device(re_cuda.cell_interps.vals.w) |> typeof <: CuDevice
  @test device(re_cuda.cell_interps.vals.ξ) |> typeof <: CuDevice
  @test device(re_cuda.cell_interps.vals.∇N_ξ) |> typeof <: CuDevice
  @test device(re_cuda.cell_interps.vals.∇∇N_ξ) |> typeof <: CuDevice
end

@testset ExtendedTestSet "AdaptExt - CPU" begin
  re = ReferenceFE(Quad4{Lagrange, 2}())
  re = Adapt.adapt_structure(Array, re)
end

@testset ExtendedTestSet "AdaptExt - CUDA" begin
  if CUDA.has_cuda()
    re = ReferenceFE(Quad4{Lagrange, 2}())
    re = Adapt.adapt_structure(CuArray, re)
    test_adapt_cuda(re)
  end
end
