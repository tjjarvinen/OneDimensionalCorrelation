
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


function reduced_ci_orbitals(b::BasisLobatto, orbitals::AbstractMatrix; pointlike=false)
    @argcheck length(b) == size(orbitals, 1)
    # Index range, for elements-2
    function _index_gen(l, i)
        if i == 0
            return 1
        else
            return l + _index_gen(l, i-1) - 1
        end
    end
    # rows for each element contributions
    function _row_range(l, i)
        return 2+(i-1)*(l-2):i*(l-2)+1 
    end
    nelements = length(b.egvector)
    new_orbitals = zeros(length(b), length(b)-nelements)
    le = length(get_element(b, 1))  # all elements have same ammount of points
    w = get_weight(b)
    if pointlike
        for ne in 1:nelements
            nr = _index_gen(le, ne-1):_index_gen(le, ne)
            tmp = diagm(ones(le))[:,begin:end-1]
            tmp[:,1] = orbitals[nr]
            t = gram_schmit(tmp, w[nr])
            new_orbitals[nr, _row_range(le, ne)] = t[:,begin+1:end]
        end
    else
        for ne in 1:nelements
            nr = _index_gen(le, ne-1):_index_gen(le, ne)
            tmp = zeros(length(b.egvector[ne]), length(b.egvector[ne])-1)
            tmp[:,1] = orbitals[nr]
            for i in 2:size(tmp,2)
                tmp[:, i] = particle_in_box(b.egvector[ne], i-1)
            end
            t = gram_schmit(tmp, w[nr])
            new_orbitals[nr, _row_range(le, ne)] = t[:,begin+1:end]
        end
    end
    new_orbitals[:,1] = orbitals[:,1]
    return gram_schmit(new_orbitals, w)
end


function boundary_ci_orbitals(b::BasisLobatto, orbitals::AbstractMatrix)
    out = zeros(length(b), number_of_elements(b)+2)
    out[:,1] = orbitals[:,1] # HF ground state
    out[1,2] = 1
    nb = 1
    for j in 1:number_of_elements(b)
        le = length( get_element(b,j) )
        nb += le - 1
        out[nb,j+2] = 1
    end
    return gram_schmit(out, get_weight(b))
end


function ci_vector_product(op, ci1, ci2; i1=1, i2=1)
    states1 = ci1["states"][:,i1]
    states2 = ci2["states"][:,i2]

    orbitals1 = ci1["orbitals"]
    orbitals2 = ci2["orbitals"]

    vector1 = ci1["indexes"]
    vector2 = ci2["indexes"]

    out = Threads.Atomic{Float64}(0)
    Threads.@threads for i in length(vector1):-1:1
        ψ1 = @view orbitals1[:, vector1[i][1] ]
        ψ2 = @view orbitals1[:, vector1[i][2] ]
        tmp = 0.0
        for j in reverse( axes(vector2,1) )
            ϕ1 = @view orbitals2[:, vector2[j][1] ]
            ϕ2 = @view orbitals2[:, vector2[j][2] ]
            tmp += op( ψ1, ψ2, ϕ1, ϕ2 ) * ( states1[i] * states2[j] )
        end
        Threads.atomic_add!(out, tmp)
    end
    return out.value
end


struct CIHamilton
    w::Vector{Float64}
    h1::Matrix{Float64}
    ve::Matrix{Float64}
    function CIHamilton(b::AbstractBasis, Vn, Ve=x->exp(-x^2))
        h1 = one_electron_operator(b, Vn)
        w = get_weight(b)

        # Electron repulsion tensor
        ve = zeros( (length(b), length(b)) )
        for j in axes(ve,2)
            ve[:,j] = Ve.(b[j].-b) 
        end
        new(w, h1, ve)
    end
end

function (cie::CIHamilton)(psi1, psi2, phi1, phi2)
    @tullio tmp = cie.w[n] * psi1[n] * phi1[n] * cie.ve[n,m] * psi2[m] * phi2[m] * cie.w[m]
    @tullio s1 = phi1[n] * psi1[n] * cie.w[n]
    @tullio s2 = phi2[n] * psi2[n] * cie.w[n]
    tmp += ( psi1' * cie.h1 * phi1 ) * s2
    tmp += ( psi2' * cie.h1 * phi2 ) * s1
    return tmp
end

struct CIOverlap
    g::Diagonal{Float64, Vector{Float64}}
    function CIOverlap(b::AbstractBasis)
        new(metric_tensor(b))
    end
end

function (cio::CIOverlap)(psi1, psi2, phi1, phi2)
    return ( psi1' * cio.g * phi1 ) * ( psi2' * cio.g * phi2 )
end