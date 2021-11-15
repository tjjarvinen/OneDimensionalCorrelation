module OneDimensionalCorrelation

using ArgCheck
using BlockBandedMatrices
using QuadGK
using LinearAlgebra
using Octavian
using PolynomialBases
using ProgressMeter
using StaticArrays
using Tullio

include("basis.jl")
include("hartree_fock.jl")
include("integrals.jl")
include("ci.jl")



export  Basis,
        BasisLobatto,
        CIHamilton,
        CIOverlap,
        CIVector,
        Element1D,
        ElementGrid,
        ElementGridLobatto

export  block_ci,
        bracket,
        ci_matrix,
        ci_vector_product,
        coulomb_matrix,
        derivative_matrix,
        element_ci_orbitals,
        eri,
        erig,
        exchange_matrix,
        fock_matrix,
        fock_matrix!,
        full_ci,
        get_element,
        get_element_indexes,
        get_identity,
        get_length,
        get_weight,
        gram_schmit,
        g_tensor,
        initial_orbitals,
        kinetic_energy,
        metric_tensor,
        one_electron_operator,
        particle_in_box,
        reduced_ci_orbitals,
        solve_hartree_fock

# Write your package code here.



end
