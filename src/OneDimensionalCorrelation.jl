module OneDimensionalCorrelation

using BlockArrays
using BlockBandedMatrices
using LinearAlgebra
using PolynomialBases

include("basis.jl")


export Element1D,
       ElementGrid,
       Basis

export derivative_matrix,
       get_weight

# Write your package code here.



end
