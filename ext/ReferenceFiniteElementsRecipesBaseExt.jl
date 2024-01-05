module ReferenceFiniteElementsRecipesBaseExt

using LaTeXStrings
using RecipesBase
using ReferenceFiniteElements
using StaticArrays

function is_inside_triangle(point)
  x_cond = (point[1] >= 0.) && (point[1] <= 1.)
  y_cond = (point[2] >= 0.) && (point[2] <= 1. - point[1])
  return x_cond && y_cond
end

@recipe function f(
  re::ReferenceFE{I, F, N, 2, Q, RefFEType, Interp, VOM, M1, M2, V}
) where {I, F, N, Q, RefFEType, Interp, VOM, M1, M2, V}
  xlabel --> L"\xi"
  ylabel --> L"\eta"
  legend --> nothing

  # draw nodes
  @series begin
    seriestype := :scatter
    color := :black
    re.nodal_coordinates[1, :], re.nodal_coordinates[2, :]
  end

  # draw edges
  edges = re.edge_nodes
  for edge in eachcol(edges)
    @series begin
      seriestype := :path
      color := :black
      temp_xs = re.nodal_coordinates[1, edge]
      temp_ys = re.nodal_coordinates[2, edge]
      temp_xs, temp_ys
    end
  end

  # draw quadrature points
  for ξ in quadrature_points(re)
    @series begin
      seriestype := :scatter
      color := :red
      [ξ[1],], [ξ[2],]
    end
  end
end

@recipe function f(
  re::ReferenceFE{I, F, N, 3, Q, RefFEType, Interp, VOM, M1, M2, V}
) where {I, F, N, Q, RefFEType, Interp, VOM, M1, M2, V}
  xlabel --> L"\xi"
  ylabel --> L"\eta"
  zlabel --> L"\zeta"
  legend --> nothing

  # draw nodes
  @series begin
    seriestype := :scatter
    color := :black
    re.nodal_coordinates[1, :], re.nodal_coordinates[2, :], re.nodal_coordinates[3, :]
  end

  # draw edges
  edges = re.edge_nodes
  for edge in eachcol(edges)
    @series begin
      seriestype := :path
      color := :black
      temp_xs = re.nodal_coordinates[1, edge]
      temp_ys = re.nodal_coordinates[2, edge]
      temp_zs = re.nodal_coordinates[3, edge]
      temp_xs, temp_ys, temp_zs
    end
  end

  # draw quadrature points
  for ξ in quadrature_points(re)
    @series begin
      seriestype := :scatter
      color := :red
      [ξ[1],], [ξ[2],], [ξ[3],]
    end
  end
end

@recipe function f(
  re::ReferenceFE{I, F, N, 2, Q, RefFEType, Interp, VOM, M1, M2, V},
  n::Int
) where {I, F, N, Q, RefFEType, Interp, VOM, M1, M2, V}
  xlabel --> L"\xi"
  ylabel --> L"\eta"
  # legend --> nothing
  # coords = re.nodal_coordinates

  # draw shape function
  if isa(re.ref_fe_type, ReferenceFiniteElements.AbstractQuad)
    ξs = -1.0:0.01:1.0
    ηs = -1.0:0.01:1.0
  elseif isa(re.ref_fe_type, ReferenceFiniteElements.AbstractTri)
    ξs = 0.0:0.01:1.0
    ηs = 0.0:0.01:1.0
  else
    @assert false
  end

  zs = Matrix{Float64}(undef, length(ξs), length(ηs))

  for i in axes(ξs, 1)
    for j in axes(ηs, 1)
      ξ = SVector{2, Float64}(ξs[i], ηs[j])
      if !is_inside_triangle(ξ)
        zs[i, j] = NaN
      else
        zs[i, j] = ReferenceFiniteElements.shape_function_values(re.ref_fe_type, SVector, ξ)[n]
      end
    end
  end

  @series begin
    seriestype := :heatmap
    color := :turbo
    colorbar := true
    ξs, ηs, zs
  end

  # draw nodes
  @series begin
    seriestype := :scatter
    color := :black
    markersize := 8
    label := nothing
    re.nodal_coordinates[1, :], re.nodal_coordinates[2, :]
  end

  # draw edges
  edges = re.edge_nodes
  for edge in eachcol(edges)
    @series begin
      seriestype := :path
      color := :black
      linewidth := 4
      label --> nothing
      temp_xs = re.nodal_coordinates[1, edge]
      temp_ys = re.nodal_coordinates[2, edge]
      temp_xs, temp_ys
    end
  end

  # draw quadrature points
  for ξ in quadrature_points(re)
    @series begin
      seriestype := :scatter
      color := :red
      markersize := 8
      label --> nothing
      [ξ[1],], [ξ[2],]
    end
  end
end

end # module
