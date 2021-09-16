# Note time to run is 30mins to 1.5 hours depending on computer

using OneDimensionalCorrelation

##
b4 = BasisLobatto(-12, 12, 4, 24)
b6 = BasisLobatto(-12, 12, 6, 16)
b8 = BasisLobatto(-12, 12, 8, 12)

# Nuclear potential
V(x) = -2( exp(-0.5(x-2)^2) + exp(-0.5(x+2)^2) )

##

E4, Eorb4, orbitals4 = solve_hartree_fock(b4, V)
E6, Eorb6, orbitals6 = solve_hartree_fock(b6, V)
E8, Eorb8, orbitals8 = solve_hartree_fock(b8, V)

##

# Takes some time ~ 3*5mins on good computer
ci4 = full_ci(b4, orbitals4, V)
ci6 = full_ci(b6, orbitals6, V)
ci8 = full_ci(b8, orbitals8, V)

## Remove boundary states

# A little faster than previous one
rorbs4 = reduced_ci_orbitals(b4, orbitals4; pointlike=true)
cir4 = full_ci(b4, rorbs4, V)

rorbs6 = reduced_ci_orbitals(b6, orbitals8; pointlike=true)
cir6 = full_ci(b6, rorbs6, V)

rorbs8 = reduced_ci_orbitals(b8, orbitals8; pointlike=true)
cir8 = full_ci(b8, rorbs8, V)

##

println("E_hf            = ", round(E6; digits=7))
println("E_fci           = ", round(ci6["energies"][1]; digits=7))

println("E_rci (4 elems) = ", round(cir4["energies"][1]; digits=7))
println("E_rci (6 elems) = ", round(cir6["energies"][1]; digits=7))
println("E_rci (8 elems) = ", round(cir8["energies"][1]; digits=7))


##

# Visualizations

using Plots

# HF ground state
plot(b6, orbitals6[:,1])

# Nuclear potential
plot(b6, V)

# Some virtual orbitals used in RCI6
# Note that orthogonalization makes them not to be
# completely point like 
plot(b6, rorbs6[:, [2,20,40,60]])