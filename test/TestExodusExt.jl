@testset ExtendedTestSet "ExodusExt" begin
  exo = ExodusDatabase("mesh_multi_block.g", "r")
  blocks = read_sets(exo, Block)
  close(exo)

  re_1 = ReferenceFE(blocks[1], 1)
  re_2 = ReferenceFE(Tet4{Lagrange, 1}())

  @test re_1.Xs == re_2.Xs
  @test re_1.face_nodes == re_2.face_nodes
  @test re_1.interior_nodes == re_2.interior_nodes
  # @test re_1.cell_interps.vals == re_2.cell_interps.vals

  @test re_1.cell_interps.vals.N == re_2.cell_interps.vals.N
  @test re_1.cell_interps.vals.∇N_ξ == re_2.cell_interps.vals.∇N_ξ
  @test re_1.cell_interps.vals.∇∇N_ξ == re_2.cell_interps.vals.∇∇N_ξ

  re_1 = ReferenceFE(blocks[2], 1)
  re_2 = ReferenceFE(Tet10{Lagrange, 1}())

  @test re_1.Xs == re_2.Xs
  @test re_1.face_nodes == re_2.face_nodes
  @test re_1.interior_nodes == re_2.interior_nodes
  # @test re_1.cell_interps.vals == re_2.cell_interps.vals
  @test re_1.cell_interps.vals.N == re_2.cell_interps.vals.N
  @test re_1.cell_interps.vals.∇N_ξ == re_2.cell_interps.vals.∇N_ξ
  @test re_1.cell_interps.vals.∇∇N_ξ == re_2.cell_interps.vals.∇∇N_ξ
end
