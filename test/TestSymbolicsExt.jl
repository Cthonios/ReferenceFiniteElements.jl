function test_symbolic_fe_fes()
    el_types = [Edge, Quad]
    # el_types = [Quad]
    for el_type in el_types
      @testset "$el_type Tests" begin
        qs = 1:10
        for q in qs
          # re = ReferenceFE{SArray, el_type, Lagrange, 1, q}()
          re = ReferenceFE{SArray, Num, el_type, Lagrange, 1, q}()
          test_q_points_inside_element(re)
          test_q_weight_positivity(re)
          test_q_weight_sum(re)
        end
  
        # just test initializing this
        # re = ReferenceFE{SArray, el_type, Lagrange, 0, 1}()
        # test_partition_of_unity_on_values(re)
        # test_partition_of_unity_on_gradients(re)
        # test_partition_of_unity_on_hessians(re)
        # test_kronecker_delta_property(re)
        
        ps = 1:5
        for p in ps
          # if p > 5
          #   arr_type = Array
          # else
          #   arr_type = SArray
          # end
          arr_type = SArray
          re = ReferenceFE{arr_type, Num, el_type, Lagrange, p, p}()
          backend = ReferenceFiniteElements.ArrayBackend{arr_type}()
          test_partition_of_unity_on_values(re)
          test_partition_of_unity_on_gradients(re)
          test_partition_of_unity_on_hessians(re)
          test_kronecker_delta_property(re, backend)
        end
      end
    end
end
