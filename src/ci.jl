
function ci_matrix(b::AbstractBasis, orbitals::AbstractMatrix, Vn; Ve=x->exp(-x^2))
    # index conversion from HF space to CI space
    indx = [ (i,j) for i in 1:length(b) for j in i:length(b) ]
    w = get_weight(b)
    ci = zeros(length(indx), length(indx))

    h1 = one_electron_operator(b, Vn)

    # Electron repulsion tensor to help calculation
    ve = zeros( eltype(orbitals), (length(b), length(b)) )
    for j in axes(ve,2)
        ve[:,j] = Ve.(b[j].-b) 
    end

    Threads.@threads for i in axes(ci,1)
        ϕ₁ = @view orbitals[:,indx[i][1]]
        ϕ₂ = @view orbitals[:,indx[i][2]]
        for j in axes(ci,2)
            ψ₁ = @view orbitals[:,indx[j][1]]
            ψ₂ = @view orbitals[:,indx[j][2]]
            # Calculate Coulomb repulsion
            @tullio tmp = w[n] * ϕ₁[n] * ψ₁[n] * ve[n,m] * ϕ₂[m] * ψ₂[m] * w[m]
            # Orbitals are orthogonal
            if indx[i][1] == indx[j][1] 
                tmp += ϕ₂' * h1 * ψ₂ 
            end
            if indx[i][2] == indx[j][2]
                tmp += ϕ₁' * h1 * ψ₁ 
            end
            ci[i,j] = tmp
        end
    end
    return ci, indx
end


function full_ci(b::AbstractBasis, orbitals::AbstractMatrix, Vn; Ve=x->exp(-x^2))
    # We need to adjust threads so that we don't oversubscribe
    bt = BLAS.get_num_threads()
    nt = Threads.nthreads()
    BLAS.set_num_threads(1)

    ci, indx = ci_matrix(b, orbitals, Vn; Ve=Ve)

    BLAS.set_num_threads(nt)

    ev, ve = eigen(ci)

    BLAS.set_num_threads(bt)

    return Dict(
        "orbitals"=>orbitals,
        "energies"=>ev,
        "states"=>ve,
        "indexes"=>indx
    )
end