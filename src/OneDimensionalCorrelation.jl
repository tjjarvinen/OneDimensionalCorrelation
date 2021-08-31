module OneDimensionalCorrelation

using ArgCheck
using BlockArrays
using BlockBandedMatrices
using QuadGK
using LinearAlgebra
using PolynomialBases

include("basis.jl")
include("wavefunction.jl")
include("integrals.jl")



export  Basis,
        BasisLobatto,
        Element1D,
        ElementGrid,
        ElementGridLobatto

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
        metric_tensor,
        particle_in_box

# Write your package code here.



end
