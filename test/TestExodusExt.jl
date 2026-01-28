function test_exodus_ext()
  exo = ExodusDatabase("mesh_multi_block.g", "r")
  blocks = read_sets(exo, Block)
  close(exo)

  re_1 = ReferenceFE(blocks[1], Lagrange; p_order = 1)
  re_2 = ReferenceFE(Tet{Lagrange, 1}(), GaussLobattoLegendre(1))

  # @test re_1.cell_Xs == re_2.cell_Xs
  # @test re_1.face_nodes == re_2.face_nodes
  # @test re_1.interior_nodes == re_2.interior_nodes
  # @test re_1.cell_interps.vals == re_2.cell_interps.vals

  @test dof_coordinates(re_1) == dof_coordinates(re_2)
  @test boundary_dofs(re_1) == boundary_dofs(re_1)
  @test interior_dofs(re_1) == interior_dofs(re_2)

  # @test re_1.cell_interps.vals.N == re_2.cell_interps.vals.N
  # @test re_1.cell_interps.vals.∇N_ξ == re_2.cell_interps.vals.∇N_ξ
  # @test re_1.cell_interps.vals.∇∇N_ξ == re_2.cell_interps.vals.∇∇N_ξ
  for q in 1:num_cell_quadrature_points(re_1)
    @test cell_shape_function_value(re_1, q) == cell_shape_function_value(re_2, q)
    @test cell_shape_function_gradient(re_1, q) == cell_shape_function_gradient(re_2, q)
    @test cell_shape_function_hessian(re_1, q) == cell_shape_function_hessian(re_2, q)
  end

  re_1 = ReferenceFE(blocks[2], Lagrange; p_order = 1)
  re_2 = ReferenceFE(Tet{Lagrange, 2}(), GaussLobattoLegendre(2))

  # @test re_1.cell_Xs == re_2.cell_Xs
  # @test re_1.face_nodes == re_2.face_nodes
  # @test re_1.interior_nodes == re_2.interior_nodes
  # @test re_1.cell_interps.vals == re_2.cell_interps.vals
  
  @test dof_coordinates(re_1) == dof_coordinates(re_2)
  @test boundary_dofs(re_1) == boundary_dofs(re_1)
  @test interior_dofs(re_1) == interior_dofs(re_2)

  # @test re_1.cell_interps.vals.N == re_2.cell_interps.vals.N
  # @test re_1.cell_interps.vals.∇N_ξ == re_2.cell_interps.vals.∇N_ξ
  # @test re_1.cell_interps.vals.∇∇N_ξ == re_2.cell_interps.vals.∇∇N_ξ
  for q in 1:num_cell_quadrature_points(re_1)
    @test cell_shape_function_value(re_1, q) == cell_shape_function_value(re_2, q)
    @test cell_shape_function_gradient(re_1, q) == cell_shape_function_gradient(re_2, q)
    @test cell_shape_function_hessian(re_1, q) == cell_shape_function_hessian(re_2, q)
  end
end
