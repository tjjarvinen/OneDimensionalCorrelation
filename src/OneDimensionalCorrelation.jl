module OneDimensionalCorrelation

using BlockArrays
using BlockBandedMatrices
using QuadGK
using LinearAlgebra
using PolynomialBases

include("basis.jl")
include("wavefunction.jl")
include("integrals.jl")



export  Element1D,
        ElementGrid,
        Basis

export  bracket,
        derivative_matrix,
        eri,
        get_identity,
        get_weight

# Write your package code here.



end
