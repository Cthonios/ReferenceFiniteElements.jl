@testset ExtendedTestSet "AdaptExt - CPU" begin
  re = ReferenceFE(Hex8(2))
  re = Adapt.adapt_structure(Array, re)

end

@testset ExtendedTestSet "AdaptExt - CUDA" begin
  if CUDA.has_cuda()
    re = ReferenceFE(Hex8(2))
    re = Adapt.adapt_structure(CuArray, re)
    test_adapt_cuda(re)
  end
end
