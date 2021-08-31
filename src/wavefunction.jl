
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