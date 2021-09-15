# OneDimensionalCorrelation

[![Build Status](https://github.com/tjjarvinen/OneDimensionalCorrelation.jl/workflows/CI/badge.svg)](https://github.com/tjjarvinen/OneDimensionalCorrelation.jl/actions)
[![Coverage](https://codecov.io/gh/tjjarvinen/OneDimensionalCorrelation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tjjarvinen/OneDimensionalCorrelation.jl)

This program is meant to study 1-dimensional electron correlation.

To install type in Julia console
```julia
using Pkg
pkg"add https://github.com/tjjarvinen/OneDimensionalCorrelation.jl"
```

Example use case
```julia
using OneDimensionalCorrelation

# Basis from x=-12 to x=12 with 4 elements and 64th order polynomials
b = BasisLobatto(-12, 12, 4, 64)

# Nuclear potential
# Two nuclears at x=-2 and x=2
V(x) = -2( exp(-0.5(x-2)^2) + exp(-0.5(x+2)^2) )

# Electron-electron potential is exp(-(x1-x2)^2)

# Solve Hartree Forck for 2 electrons
E, Eorb, orbitals = solve_hartree_fock(b, V)

# Solve Full CI
ci = full_ci(b, orbitals, V)
```

To visualize results
```julia
using Plots

# Plot nuclear potential
plot(b, V)

# Plot 4 lowest orbitals
plot(b, orbitals[:,1:4])
```