
function bracket(
        b::Basis,
        psi1::AbstractVector,
        op::AbstractMatrix,
        psi2::AbstractVector
    )
    w = get_weight(b)
    return (conj.(psi1).*w)' * op * psi2
end

function bracket(b::Basis, psi1::AbstractVector, psi2::AbstractVector)
    return bracket(b, psi1, get_identity(b), psi2)
end