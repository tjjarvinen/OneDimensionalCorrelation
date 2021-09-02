
"""
    erig(b::Basis, i::Int, j::Int; alpha=1.0, scale=1.0)
    erig(b::Basis, i::Int, j::Int, n::Int, m::Int; alpha=1.0, scale=1.0)

Electron repulsion integral with Gaussian repulsion scale*exp(-α*r^2).

For erig4 `i` and `j` are indices for particale 1 and `n` and `m` for particale 2.
"""
function erig(b::AbstractBasis, i::Int, j::Int; Ve=x->exp(-x^2))
    return Ve(b[i]-b[j])
end


function erig(b::AbstractBasis, i::Int, j::Int, n::Int, m::Int; Ve=x->exp(-x^2))::Float64
    # TODO this function might not be type stable
    T = Float64   # This need to change for AD to work eg. Base.return_type()
    # Basis functions are orthogonal
    i != j && return zero(T)
    n != m && return zero(T)
    return Ve(b[i]-b[n])
end


function fock_matrix(b::AbstractBasis; Ve=x->exp(-x^2), Vn=x->-2exp(-0.5x^2))
    orbitals = initial_orbitals(b)
    return fock_matrix(b, orbitals; Ve=Ve, Vn=Vn)
end


function fock_matrix(b::AbstractBasis, orbitals::AbstractMatrix; Ve=x->exp(-x^2), Vn=x->-2exp(-0.5x^2))
    @argcheck size(orbitals) == (length(b), length(b))
    h₁ = one_electron_operator(b, Vn)
    J = coulomb_matrix(b, orbitals)
    #K = exchange_matrix(b, orbitals) # For more than 2 electrons
    return h₁ + J  # h₁ + 2*J - K
end

function fock_matrix!(f::AbstractMatrix, b::AbstractBasis, orbitals::AbstractMatrix; Ve=x->exp(-x^2), Vn=x->-2exp(-0.5x^2))
    @argcheck size(f) == size(orbitals) == (length(b), length(b))
    h₁ = one_electron_operator(b, Vn)
    J = coulomb_matrix(b, orbitals)
    #K = exchange_matrix(b, orbitals) # For more than 2 electrons
    return f .= h₁ + J  # h₁ + 2*J - K
end


function coulomb_matrix(b::AbstractBasis, orbitals::AbstractMatrix)
    C = zeros(size(orbitals))
    return coulomb_matrix!(C, b, orbitals)
end

function coulomb_matrix!(C::AbstractMatrix, b::AbstractBasis, orbitals::AbstractMatrix)
    @argcheck size(C) == size(orbitals) == (length(b), length(b))
    l = length(b)
    w = get_weight(b)

    # Two electrons in total
    ρ = orbitals[:,1] * orbitals[:,1]'
    for i in 1:l
        # C[i,j]=0 for i!=j
        # also erig(b,i,i,n,m) = 0 for n!=m
        C[i,i] = sum( n -> erig(b,i,i,n,n) * ρ[n,n] * w[n], 1:l)
        C[i,i] *= w[i]
    end
    return C
end


function exchange_matrix(b::AbstractBasis, orbitals::AbstractMatrix)
    @argcheck size(orbitals) == (length(b), length(b))
    l = length(b)
    K = zeros(l,l)
    w = get_weight(b)

    # Two electrons in total
    ρ = orbitals[:,1] * orbitals[:,1]'
    Threads.@threads for i in 1:l
        for j in 1:l
            for n in 1:l
                K[i,j] += sum( k-> erig(b, i,n,j,k) * ρ[n,k], 1:l)
            end
            K[i,j] *= w[i] * w[j]
        end
    end
    return K     
end


function bracket(
        b::AbstractBasis,
        psi1::AbstractVector,
        op::AbstractMatrix,
        psi2::AbstractVector
    )
    w = get_weight(b)
    return (conj.(psi1).*w)' * op * psi2
end

function bracket(b::AbstractBasis, psi1::AbstractVector, psi2::AbstractVector)
    w = get_weight(b)
    return sum( conj.(psi1) .* w .* psi2 )
end


function kinetic_energy(b::AbstractBasis; mass=1.0)
    ∇ = derivative_matrix(b)
    g = metric_tensor(b)
    return (0.5/mass) * ∇' * g * ∇ # = 0.5 * (g∇)' * ∇
end


function one_electron_operator(b, Vn; mass=1.0)
    g = metric_tensor(b)
    V = g * diagm( Vn.(b) )
    return kinetic_energy(b; mass=mass) + V
end