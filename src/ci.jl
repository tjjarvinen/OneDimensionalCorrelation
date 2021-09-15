
function ci_matrix(b::AbstractBasis, orbitals::AbstractMatrix, Vn; Ve=x->exp(-x^2))
    @argcheck length(b) == size(orbitals, 1)
    # index conversion from HF space to CI space
    indx = [ (i,j) for i in axes(orbitals, 2) for j in axes(orbitals, 2) ]
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

    @info "Creating CI matrix"
    @time ci, indx = ci_matrix(b, orbitals, Vn; Ve=Ve)

    BLAS.set_num_threads(nt)

    @info "Diagonalizing CI matrix"
    @time ev, ve = eigen(ci)

    BLAS.set_num_threads(bt)

    return Dict(
        "orbitals"=>orbitals,
        "energies"=>ev,
        "states"=>ve,
        "indexes"=>indx
    )
end


function reduced_ci_orbitals(b::BasisLobatto, orbitals::AbstractMatrix)
    @argcheck length(b) == size(orbitals, 1)
    function _gram_schmit(orbitals::AbstractMatrix, w::AbstractVector)
        # This is a special version of Gram-Schmid
        # Where 1st orbital is ignored from normalization
        # This is becaus it is expected to 
        out = similar(orbitals)
        norm = sqrt( sum( w .* orbitals[:,1].^2 ) )
        out[:,1] = orbitals[:,1] # 1st on is expected to be normalized
        for i in 2:size(out,2)
            tmp = orbitals[:,i]
            for j in 1:i-1
               tmp -= sum( out[:,j] .* w .* orbitals[:,i] ) .* out[:,j]   
            end
            norm = sqrt( sum( w .* tmp.^2 ) )
            out[:,i] = tmp ./ norm 
        end
        return out
    end
    function _index_gen(l, i)
        if i == 0
            return 1
        else
            return l + _index_gen(l, i-1) - 1
        end
    end
    function _row_range(l, i)
        return 2+(i-1)*(l-2):i*(l-2)+1 
    end
    nelements = length(b.egvector)
    new_orbitals = zeros(length(b), length(b)-nelements)
    le = length(get_element(b, 1))  # all elements have same ammount of points
    w = get_weight(b)
    for ne in 1:nelements
        nr = _index_gen(le, ne-1):_index_gen(le, ne)
        tmp = diagm(ones(le))[:,begin:end-1]
        tmp[:,1] = orbitals[nr]
        t = _gram_schmit(tmp, w[nr])
        new_orbitals[nr, _row_range(le, ne)] = t[:,begin+1:end]
    end
    new_orbitals[:,1] = orbitals[:,1]
    return new_orbitals
end

