module OneDimensionalCorrelation

using ArgCheck
using QuadGK
using LinearAlgebra
using PolynomialBases
using ProgressMeter
using Tullio

include("basis.jl")
include("hartree_fock.jl")
include("integrals.jl")
include("ci.jl")



export  Basis,
        BasisLobatto,
        Element1D,
        ElementGrid,
        ElementGridLobatto

export  bracket,
        ci_matrix,
        coulomb_matrix,
        derivative_matrix,
        eri,
        erig,
        exchange_matrix,
        fock_matrix,
        fock_matrix!,
        full_ci,
        get_element,
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
