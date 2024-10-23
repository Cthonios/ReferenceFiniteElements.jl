module ReferenceFiniteElementsRecipesBaseExt

using LaTeXStrings
using RecipesBase
using ReferenceFiniteElements
using StaticArrays

@recipe function f(re::ReferenceFE{I, F, <:Edge}, n::Int) where {I, F}
  xlabel --> L"\xi"
  ylabel --> L"\eta"
  legend --> nothing

  # draw nodes
  @series begin
    seriestype := :scatter
    color := :black
    temp = vcat(map(x -> vcat(x, 0.), re.Xs)'...)
    temp[:, 1], temp[:, 2]
  end

  @series begin
    seriestype := :scatter
    color := :red
    temp = vcat(map(x -> vcat(x, 0.), re.ξs)'...)
    temp[:, 1], temp[:, 2]
  end

  ξs = LinRange(-1., 1., 100)
  backend = ReferenceFiniteElements.ArrayBackend{SArray}()
  Ns = shape_function_value.((re.element,), (re.Xs,), ξs, (backend,))
  Ns = map(x -> x[n], Ns)

  @series begin
    seriestype := :path
    color := :blue
    ξs, Ns
  end
end

@recipe function f(re::ReferenceFE{I, F, <:Edge}, n::Int, d::Int) where {I, F}
  xlabel --> L"\xi"
  ylabel --> L"\eta"
  legend --> nothing

  # draw nodes
  @series begin
    seriestype := :scatter
    color := :black
    temp = vcat(map(x -> vcat(x, 0.), re.Xs)'...)
    temp[:, 1], temp[:, 2]
  end

  @series begin
    seriestype := :scatter
    color := :red
    temp = vcat(map(x -> vcat(x, 0.), re.ξs)'...)
    temp[:, 1], temp[:, 2]
  end

  ξs = LinRange(-1., 1., 100)
  backend = ReferenceFiniteElements.ArrayBackend{SArray}()
  Ns = shape_function_gradient.((re.element,), (re.Xs,), ξs, (backend,))
  Ns = map(x -> x[n], Ns)

  display(Ns)
  @series begin
    seriestype := :path
    color := :blue
    ξs, Ns
  end
end

end # module