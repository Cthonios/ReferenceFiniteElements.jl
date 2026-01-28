module RecipesBaseExt

using RecipesBase
using ReferenceFiniteElements

@recipe function f(re::ReferenceFE, plot_type, index_1, index_2)
  if dimension(re) == 1
    xlabel --> "x"
    if plot_type == :shape_function_value
      ylabel --> "N"
    elseif plot_type == :shape_function_gradient
      ylabel --> "∇N_ξ"
    else
      @assert false
    end

    # draw vertices
    @series begin 
      seriestype := :scatter
      color := :black
      temp = vertex_coordinates(re)
      temp[1, :], temp[2, :]
    end

    Xs = dof_coordinates(re)
    @series begin
      seriestype := :scatter
      color := :blue
      Xs[1, :], zeros(length(Xs[1, :]))
    end
  
    # evaluate values at certain points
    if re.element.shifted
      ξs = LinRange(0., 1., 101)
    else
      ξs = LinRange(-1., 1., 101)
    end

    if plot_type == :shape_function_value
      Ns = mapreduce(x -> ReferenceFiniteElements.shape_function_value(re.element, Xs, x), hcat, ξs)
    elseif plot_type == :shape_function_gradient
      Ns = map(x -> ReferenceFiniteElements.shape_function_gradient(re.element, Xs, x), ξs)
      Ns = mapreduce(x -> x[index_2, :], hcat, Ns)
    else
      @assert false
    end

    for n in axes(Ns, 1)
      @series begin
        seriestype := :path
        ξs, Ns[n, :]
      end
    end
  else
    @assert false "Unsupported dimension for plotting $(dimension(re))"
  end
end

# @recipe function f(re::ReferenceFE{I, F, <:Edge}, n::Int) where {I, F}
#   xlabel --> L"\xi"
#   ylabel --> L"\eta"
#   legend --> nothing

#   # draw nodes
#   @series begin
#     seriestype := :scatter
#     color := :black
#     temp = vcat(map(x -> vcat(x, 0.), re.Xs)'...)
#     temp[:, 1], temp[:, 2]
#   end

#   @series begin
#     seriestype := :scatter
#     color := :red
#     temp = vcat(map(x -> vcat(x, 0.), re.ξs)'...)
#     temp[:, 1], temp[:, 2]
#   end

#   ξs = LinRange(-1., 1., 100)
#   backend = ReferenceFiniteElements.ArrayBackend{SArray}()
#   Ns = shape_function_value.((re.element,), (re.Xs,), ξs, (backend,))
#   Ns = map(x -> x[n], Ns)

#   @series begin
#     seriestype := :path
#     color := :blue
#     ξs, Ns
#   end
# end

# @recipe function f(re::ReferenceFE{I, F, <:Edge}, n::Int, d::Int) where {I, F}
#   xlabel --> L"\xi"
#   ylabel --> L"\eta"
#   legend --> nothing

#   # draw nodes
#   @series begin
#     seriestype := :scatter
#     color := :black
#     temp = vcat(map(x -> vcat(x, 0.), re.Xs)'...)
#     temp[:, 1], temp[:, 2]
#   end

#   @series begin
#     seriestype := :scatter
#     color := :red
#     temp = vcat(map(x -> vcat(x, 0.), re.ξs)'...)
#     temp[:, 1], temp[:, 2]
#   end

#   ξs = LinRange(-1., 1., 100)
#   backend = ReferenceFiniteElements.ArrayBackend{SArray}()
#   Ns = shape_function_gradient.((re.element,), (re.Xs,), ξs, (backend,))
#   Ns = map(x -> x[n], Ns)

#   display(Ns)
#   @series begin
#     seriestype := :path
#     color := :blue
#     ξs, Ns
#   end
# end

end # module