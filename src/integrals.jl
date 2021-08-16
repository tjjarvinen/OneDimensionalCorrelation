
function g_tensor(b::Basis)
    #TODO
end

function eri(b::Basis, i::Int, j::Int, n::Int, m::Int; tmax=1000)
    # Basis functions are orthogonal
    i != j && return 0
    n != m && return 0
    return eri(b, i, n; tmax=tmax)
end

function eri(b::Basis, i::Int, j::Int; tmax=1000, threshold=0.01)
    f(t, x) = 2/√π *exp(-t^2*x^2)
    r = abs(b[i]-b[j])
    r > threshold && return 1/r
    return quadgk(t->f(t,r), 0, tmax; rtol=1e-12)[1]
end