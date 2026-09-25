# Element types
Below is the currently implemented set of elements organized
by vertices, edges, faces, and volumes respectively. If you would
like to see an additional element type supported, please file and issue
or open a PR.

## 0-Dimensional Elements
```@docs
Vertex
```

## 1-Dimensional Elements
```@docs
Edge
```

## 2-Dimensional Elements
```@docs
Quad
Tri
```

## Enriched elements
`Tri{EnrichedLagrange, 2}` and `Tet{EnrichedLagrange, 2}` carry the degree-2
Lagrange space enriched with bubbles in a nodal basis; see `EnrichedLagrange`.

## 3-Dimensional Elements
```@docs
Hex
Tet
```

# Polynomial types
```@docs
EnrichedLagrange
Hermite
Lagrange
ReferenceFiniteElements.NedelecFirstKind
ReferenceFiniteElements.NedelecSecondKind
ReferenceFiniteElements.NoInterpolation
RaviartThomas
Serendipity
```

# Quadrature types
```@docs
GaussLegendre
GaussLobattoLegendre
```