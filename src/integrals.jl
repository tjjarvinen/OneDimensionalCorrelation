
"""
    erig(b::Basis, i::Int, j::Int; alpha=1.0, scale=1.0)
    erig(b::Basis, i::Int, j::Int, n::Int, m::Int; alpha=1.0, scale=1.0)

Electron repulsion integral with Gaussian repulsion scale*exp(-α*r^2).

For erig4 `i` and `j` are indices for particale 1 and `n` and `m` for particale 2.
"""
function erig(b::AbstractBasis, i::Int, j::Int; alpha=1.0, scale=1.0)
    return scale * exp( -alpha * (b[i]-b[j])^2 )
end


function erig(b::AbstractBasis, i::Int, j::Int, n::Int, m::Int; alpha=1.0, scale=1.0)
    # TODO this function might not be type stable
    # Basis functions are orthogonal
    T = promote_type(typeof(alpha), typeof(scale))
    i != j && return zero(T)
    n != m && return zero(T)
    return erig(b, i, n; alpha=alpha, scale=scale)
end


function fock_matrix(b::AbstractBasis; alpha=1.0, scale=1.0)
    orbitals = initial_orbitals(b)
    return fock_matrix(b, orbitals; alpha=alpha, scale=scale)
end


function fock_matrix(b::AbstractBasis, orbitals::AbstractMatrix; alpha=1.0, scale=1.0)
    ∇ = derivative_matrix(b)
    g = metric_tensor(b)
    Eₖ = 0.5 * ∇' * g * ∇
    # Two electrons so nuclear potential is 2
    Vₙₑ = -2 .* scale .*  g * diagm( exp.( -alpha .* b.^2 ) )
    J = coulomb_matrix(b, orbitals)
    K = exchange_matrix(b, orbitals)
    h₁ = Eₖ +  Vₙₑ
    return h₁ + J - K
end


function coulomb_matrix(b::AbstractBasis, orbitals::AbstractMatrix)
    l = length(b)
    C = zeros(l, l)
    w = get_weight(b)

    # Two electrons in total
    ρ = orbitals[:,1] * orbitals[:,1]'
    for i in 1:l
        for j in 1:l
            C[i,j] = sum( n -> erig(b, i, j, n, n) * ρ[n,n] * w[n], 1:l)
            C[i,j] *= w[i]
        end
    end
    return C
end


function exchange_matrix(b::AbstractBasis, orbitals::AbstractMatrix)
    l = length(b)
    K = zeros(l,l)
    w = get_weight(b)

    # Two electrons in total
    ρ = orbitals[:,1] * orbitals[:,1]'
    for i in 1:l
        for j in 1:l
            for n in 1:l, m in 1:l
                K[i,j] += erig(b, i,n,j,m) * ρ[n,m]
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