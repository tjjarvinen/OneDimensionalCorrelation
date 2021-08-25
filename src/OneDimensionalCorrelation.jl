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
        coulomb_matrix,
        derivative_matrix,
        eri,
        erig,
        exchange_matrix,
        fock_matrix,
        get_identity,
        get_length,
        get_weight,
        g_tensor,
        initial_orbitals,
        particle_in_box

# Write your package code here.



end
