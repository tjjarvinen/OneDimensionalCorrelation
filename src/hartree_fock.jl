
function particle_in_box(b::AbstractBasis, n::Int)
    L = get_length(b)
    N = sqrt(2/L)
    a = n * π / L
    return N .* sin.(a .* (b .+ 0.5*L) )
end

function initial_orbitals(b::AbstractBasis)
    l = length(b)
    ψ = zeros(l,l)
    for j in 1:l
        ψ[:,j] = particle_in_box(b, j)
    end
    return ψ
end

function solve_hartree_fock(b::AbstractBasis, Vn ; Ve=x->exp(-x^2), rtol=1E-9, max_iter=100)
    orbitals = initial_orbitals(b)
    F = fock_matrix(b, orbitals; Vn=Vn)
    h₁ = one_electron_operator(b, Vn)
    # Energy is h₁ + F (=2h₁+J) for two electrons
    E₀ = orbitals[:,1]' * ( h₁ + F ) * orbitals[:,1]
    @info "Energy at start = $(E₀)"
    eo, orbitals = eigen(F)
    for i in 1:max_iter
        F = fock_matrix(b, orbitals; Vn=Vn)
        E₁ = orbitals[:,1]' * ( h₁ + F ) * orbitals[:,1]
        @info "Energy at iteration $i = $(E₁)"
        if abs(E₁-E₀) < rtol
            break
        else
            E₀ = E₁ 
            eo, orbitals = eigen(F)
        end
    end
    return eo, orbitals
end