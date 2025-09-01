@testset ExtendedTestSet "ExodusExt" begin
  exo = ExodusDatabase("mesh_multi_block.g", "r")
  blocks = read_sets(exo, Block)
  close(exo)

  re_1 = ReferenceFE(blocks[1], 1)
  re_2 = ReferenceFE(Tet4{Lagrange, 1}())

  @test re_1.cell_Xs == re_2.cell_Xs
  @test re_1.face_nodes == re_2.face_nodes
  @test re_1.interior_nodes == re_2.interior_nodes
  # @test re_1.cell_interps.vals == re_2.cell_interps.vals

  # @test re_1.cell_interps.vals.N == re_2.cell_interps.vals.N
  # @test re_1.cell_interps.vals.∇N_ξ == re_2.cell_interps.vals.∇N_ξ
  # @test re_1.cell_interps.vals.∇∇N_ξ == re_2.cell_interps.vals.∇∇N_ξ
  for q in 1:num_quadrature_points(re_1)
    @test shape_function_value(re_1, q) == shape_function_value(re_2, q)
    @test shape_function_gradient(re_1, q) == shape_function_gradient(re_2, q)
    @test shape_function_hessian(re_1, q) == shape_function_hessian(re_2, q)
  end

  re_1 = ReferenceFE(blocks[2], 1)
  re_2 = ReferenceFE(Tet10{Lagrange, 1}())

  @test re_1.cell_Xs == re_2.cell_Xs
  @test re_1.face_nodes == re_2.face_nodes
  @test re_1.interior_nodes == re_2.interior_nodes
  # @test re_1.cell_interps.vals == re_2.cell_interps.vals
  
  # @test re_1.cell_interps.vals.N == re_2.cell_interps.vals.N
  # @test re_1.cell_interps.vals.∇N_ξ == re_2.cell_interps.vals.∇N_ξ
  # @test re_1.cell_interps.vals.∇∇N_ξ == re_2.cell_interps.vals.∇∇N_ξ
  for q in 1:num_quadrature_points(re_1)
    @test shape_function_value(re_1, q) == shape_function_value(re_2, q)
    @test shape_function_gradient(re_1, q) == shape_function_gradient(re_2, q)
    @test shape_function_hessian(re_1, q) == shape_function_hessian(re_2, q)
  end
end
