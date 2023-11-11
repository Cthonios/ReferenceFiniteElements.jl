@testset ExtendedTestSet "ExodusExt" begin
  exo = ExodusDatabase("mesh_multi_block.g", "r")
  blocks = read_sets(exo, Block)
  close(exo)

  re_1 = ReferenceFE(blocks[1], 1)
  re_2 = ReferenceFE(Tet4(1))

  @test re_1.nodal_coordinates == re_2.nodal_coordinates
  @test re_1.face_nodes == re_2.face_nodes
  @test re_1.interior_nodes == re_2.interior_nodes
  @test re_1.interpolants == re_2.interpolants

  re_1 = ReferenceFE(blocks[2], 1)
  re_2 = ReferenceFE(Tet10(1))

  @test re_1.nodal_coordinates == re_2.nodal_coordinates
  @test re_1.face_nodes == re_2.face_nodes
  @test re_1.interior_nodes == re_2.interior_nodes
  @test re_1.interpolants == re_2.interpolants
end