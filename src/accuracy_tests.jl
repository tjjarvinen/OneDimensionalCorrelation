"""
    ref_int(a=1, b=1; tmax=10000, rtol=1E-12)

Calculates Coulomb repulsion between two Gaussian charge distributions
that have same origin.

# Arguments
- `a=1`  :   with of first Gaussian exp(-ax^2)
- `b=1`  :   with of second Gaussian exp(-bx^2)

# Keywords
- `tmax=10000`   :  maximum t-value for integration
- `rtol=1E-12`   :  relative tolerance for t-integration
"""
function ref_int(a=1, b=1; tmax=10000, rtol=1E-12)
    _f(t) = 1/sqrt( (a+b)*t^2 + a*b)
    norm = a*b/π
    return 2*√π * norm * quadgk(_f, 0, tmax; rtol=rtol)[1]
end


function test_eri_accuracy(d, ne, np; a=1, b=1)
    basis = Basis(-0.5*d, 0.5*d, ne, np)
    ρ1 = exp.(-a.*basis.^2) .* sqrt(a/π)
    ρ2 = exp.(-b.*basis.^2) .* sqrt(b/π)
    ω = get_weight(basis)
    w = diagm(ω)
    l = length(basis)

    J = zeros(l,l)
    for i in 1:l, j in 1:l
        J[i,j] = eri(basis, i, j)
    end
    cal = ρ2'*w*J*ρ1

    ref = ref_int(a, b)

    return Dict("ref"=>ref, "cal"=>cal)
end