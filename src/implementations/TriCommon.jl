abstract type AbstractTri <: ReferenceFE end

"""
Eventually move this to FastQuassQuadrature implementation
"""
function Quadrature(::E, degree::I, Rtype::Type = Float64) where {I <: Integer, E <: AbstractTri}
  if degree == 1
    # this is dumb that I have to do this to minimize allocations
    ξ = Matrix{Rtype}(undef, 2, 1)
    ξ[1, 1] = 1. / 3.
    ξ[2, 1] = 1. / 3.
    w = [0.5]
  elseif degree == 2
    ξ = [
      2. / 3. 1. / 6. 1. / 6.;
      1. / 6. 2. / 3. 1. / 6
    ]
    w = [1. / 6., 1. / 6., 1. / 6.]
  elseif degree <= 4
    ξ = [
      1.081030181680700E-01 4.459484909159650E-01 4.459484909159650E-01 8.168475729804590E-01 9.157621350977100E-02 9.157621350977100E-02;
      4.459484909159650E-01 1.081030181680700E-01 4.459484909159650E-01 9.157621350977100E-02 8.168475729804590E-01 9.157621350977100E-02
    ]

    w = [
      1.116907948390055E-01,
      1.116907948390055E-01,
      1.116907948390055E-01,
      5.497587182766100E-02,
      5.497587182766100E-02,
      5.497587182766100E-02
    ]
  elseif degree <= 5
    ξ = [
      3.33333333333333E-01 5.97158717897700E-02 4.70142064105115E-01 4.70142064105115E-01 7.97426985353087E-01 1.01286507323456E-01 1.01286507323456E-01;
      3.33333333333333E-01 4.70142064105115E-01 5.97158717897700E-02 4.70142064105115E-01 1.01286507323456E-01 7.97426985353087E-01 1.01286507323456E-01
    ]

    w = [
      1.12500000000000E-01,
      6.61970763942530E-02,
      6.61970763942530E-02,
      6.61970763942530E-02,
      6.29695902724135E-02,
      6.29695902724135E-02,
      6.29695902724135E-02
    ]
  elseif degree <= 6
    ξ = [
      5.01426509658179E-01 2.49286745170910E-01 2.49286745170910E-01 8.73821971016996E-01 6.30890144915020E-02 6.30890144915020E-02 5.31450498448170E-02 6.36502499121399E-01 3.10352451033784E-01 5.31450498448170E-02 6.36502499121399E-01 3.10352451033784E-01;
      2.49286745170910E-01 5.01426509658179E-01 2.49286745170910E-01 6.30890144915020E-02 8.73821971016996E-01 6.30890144915020E-02 3.10352451033784E-01 5.31450498448170E-02 6.36502499121399E-01 6.36502499121399E-01 3.10352451033784E-01 5.31450498448170E-02
    ]

    w = [
      5.83931378631895E-02,
      5.83931378631895E-02,
      5.83931378631895E-02,
      2.54224531851035E-02,
      2.54224531851035E-02,
      2.54224531851035E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02,
      4.14255378091870E-02
    ]
  else
    throw(ErrorException("Unsupported quadrature degree."))
  end
  return Quadrature{Rtype}(ξ, w)
end