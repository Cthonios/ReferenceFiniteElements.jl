@testset ExtendedTestSet "SimplexTri" begin

  # e = SimplexTri(1)
  # coords, face, ints = ReferenceFiniteElements.element_stencil(e, Int64, Float64)
  

  for int_type in [Int32, Int64]
    for float_type in [Float32, Float64]
      for array_type in [SArray, MArray]

        # degree 1
        e = ReferenceFE(SimplexTri(1); int_type=int_type, float_type=float_type, array_type=array_type)

        @test e.nodal_coordinates[1, 1] ≈ 1.
        @test e.nodal_coordinates[2, 1] ≈ 0.
        @test e.nodal_coordinates[1, 2] ≈ 0.
        @test e.nodal_coordinates[2, 2] ≈ 1.
        @test e.nodal_coordinates[1, 3] ≈ 0.
        @test e.nodal_coordinates[2, 3] ≈ 0.

        @test e.interpolants.N[1] ≈ [1. / 3., 1. / 3., 1. / 3.]
        @test e.interpolants.∇N_ξ[1] ≈ [1. 0.; 0. 1.; -1. -1.]

        # degree 2
        e = ReferenceFE(SimplexTri(2); int_type=int_type, float_type=float_type, array_type=array_type)
        @test e.nodal_coordinates[1, 1] ≈ 1.
        @test e.nodal_coordinates[2, 1] ≈ 0.
        @test e.nodal_coordinates[1, 3] ≈ 0.
        @test e.nodal_coordinates[2, 3] ≈ 1.
        @test e.nodal_coordinates[1, 6] ≈ 0.
        @test e.nodal_coordinates[2, 6] ≈ 0.
        #
        # @test e.nodal_coordinates[1, 2] ≈ 

        @test e.interpolants.N[1] ≈ [0.22222222, 0.44444444, -0.11111111, 0.44444444, 0.11111111, -0.11111111]
        @test e.interpolants.N[2] ≈ [-0.11111111, 0.44444444, 0.22222222, 0.11111111, 0.44444444, -0.11111111]
        @test e.interpolants.N[3] ≈ [-0.11111111, 0.11111111, -0.11111111, 0.44444444, 0.44444444, 0.22222222]

        if float_type == Float32
          tol = 1e-6
        else
          tol = 1e-14
        end

        @test e.interpolants.∇N_ξ[1][1, 1] ≈ 1.66666667e+00
        @test e.interpolants.∇N_ξ[1][1, 2] ≈ 0.0 atol=tol
        @test e.interpolants.∇N_ξ[1][2, 1] ≈ 6.66666667e-01
        @test e.interpolants.∇N_ξ[1][2, 2] ≈ 2.66666667e+00
        @test e.interpolants.∇N_ξ[1][3, 1] ≈ 0.0 atol=tol
        @test e.interpolants.∇N_ξ[1][3, 2] ≈ -1. / 3.
        @test e.interpolants.∇N_ξ[1][4, 1] ≈ -2.
        @test e.interpolants.∇N_ξ[1][4, 2] ≈ -2.66666667e+00
        @test e.interpolants.∇N_ξ[1][5, 1] ≈ -2. / 3.
        @test e.interpolants.∇N_ξ[1][5, 2] ≈ 0. atol=tol
        @test e.interpolants.∇N_ξ[1][6, 1] ≈ 1. / 3.
        @test e.interpolants.∇N_ξ[1][6, 2] ≈ 1. / 3.
        #
        @test e.interpolants.∇N_ξ[2][1, 1] ≈ -1. / 3.
        @test e.interpolants.∇N_ξ[2][1, 2] ≈ 0.0 atol=tol
        @test e.interpolants.∇N_ξ[2][2, 1] ≈ 2.66666667e+00
        @test e.interpolants.∇N_ξ[2][2, 2] ≈ 2. / 3.
        @test e.interpolants.∇N_ξ[2][3, 1] ≈ 0.0 atol=tol
        @test e.interpolants.∇N_ξ[2][3, 2] ≈ 5. / 3.
        @test e.interpolants.∇N_ξ[2][4, 1] ≈ 0. atol=tol
        @test e.interpolants.∇N_ξ[2][4, 2] ≈ -2. / 3.
        @test e.interpolants.∇N_ξ[2][5, 1] ≈ -8. / 3.
        @test e.interpolants.∇N_ξ[2][5, 2] ≈ -2.
        @test e.interpolants.∇N_ξ[2][6, 1] ≈ 1. / 3.
        @test e.interpolants.∇N_ξ[2][6, 2] ≈ 1. / 3.
        #
        @test e.interpolants.∇N_ξ[3][1, 1] ≈ -1. / 3.
        @test e.interpolants.∇N_ξ[3][1, 2] ≈ 0.0 atol=tol
        @test e.interpolants.∇N_ξ[3][2, 1] ≈ 2. / 3.
        @test e.interpolants.∇N_ξ[3][2, 2] ≈ 2. / 3.
        @test e.interpolants.∇N_ξ[3][3, 1] ≈ 0.0 atol=tol
        @test e.interpolants.∇N_ξ[3][3, 2] ≈ -1. / 3.
        @test e.interpolants.∇N_ξ[3][4, 1] ≈ 2.
        @test e.interpolants.∇N_ξ[3][4, 2] ≈ -2. / 3.
        @test e.interpolants.∇N_ξ[3][5, 1] ≈ -2. / 3.
        @test e.interpolants.∇N_ξ[3][5, 2] ≈ 2.
        @test e.interpolants.∇N_ξ[3][6, 1] ≈ -5. / 3.
        @test e.interpolants.∇N_ξ[3][6, 2] ≈ -5. / 3.
      end
    end
  end
end