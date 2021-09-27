using OneDimensionalCorrelation
##
println("Calculate CI by fist calculating CI on smaller element blocks")
println("and then combining lowest 2 CI vectors to form a new CI calculation.")
println("The result is then compared to full CI calculation.")
##
b4 = BasisLobatto(-12, 12, 4, 24)
b6 = BasisLobatto(-12, 12, 6, 16)
b8 = BasisLobatto(-12, 12, 8, 12)

# Nuclear potential
V(x) = -2( exp(-0.5(x-2)^2) + exp(-0.5(x+2)^2) )

##
println("\n *** Hartree Fock calculations ***\n")

E4, Eorb4, orbitals4 = solve_hartree_fock(b4, V)
E6, Eorb6, orbitals6 = solve_hartree_fock(b6, V)
E8, Eorb8, orbitals8 = solve_hartree_fock(b8, V)

##

println("\n *** CI in blocks calculations ***\n")

eb4 = block_ci(b4, orbitals4, V; nstates=2)
eb6 = block_ci(b6, orbitals6, V; nstates=2)
eb8 = block_ci(b8, orbitals8, V; nstates=2)


##
println("\n *** Full CI calculation ***\n")
ci = full_ci(b6, orbitals6, V)

##

E_corr = ci["energies"][1] - E6

Eb4_corr = eb4["energy"][1] - E4
Eb6_corr = eb6["energy"][1] - E6
Eb8_corr = eb8["energy"][1] - E8


E2 = ci["energies"][2]

E2_err4 = abs( ( eb4["energy"][2] - E2 ) / E2 )
E2_err6 = abs( ( eb6["energy"][2] - E2 ) / E2 )
E2_err8 = abs( ( eb8["energy"][2] - E2 ) / E2 )

##
println("\nHartree-Fock energy: ", E6)
println("Total correlation energy (ground state): ", E_corr, "\n")

println("Relative correlation energy based on number of elements when calculating with blocks:\n")
println("4 elements ", round(100 * (Eb4_corr/E_corr); digits=2),"% of total correlation")
println("6 elements ", round(100 * (Eb6_corr/E_corr); digits=2),"% of total correlation")
println("8 elements ", round(100 * (Eb8_corr/E_corr); digits=2),"% of total correlation")

println("\n***\n")

println("First excited state energy FCI: ", E2,"\n")
println("Excited state energies when using element blocks:")
println("4 elements: ", eb4["energy"][2], "  relative accuracy: ", round(E2_err4; sigdigits=2))
println("6 elements: ", eb6["energy"][2], "  relative accuracy: ", round(E2_err6; sigdigits=2))
println("8 elements: ", eb8["energy"][2], "  relative accuracy: ", round(E2_err8; sigdigits=2))