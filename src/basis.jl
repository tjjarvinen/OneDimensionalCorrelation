

abstract type AbstractElementGrid{T} <: AbstractVector{T} end
abstract type AbstractBasis{T} <: AbstractVector{T} end

struct Element1D
    a::Float64
    b::Float64
    function Element1D(a,b)
        @assert b>a "Element1D creation error, b not >a"
        new(a,b)
    end
end

struct ElementGrid <: AbstractElementGrid{Float64}
    el::Element1D
    x::Vector{Float64}
    w::Vector{Float64}
    function ElementGrid(a, b, n)
        x, w = gausslegendre(n)
        s = (a+b)
        xx = (b-a)/2 .* x .+ (a+b)/2
        new(Element1D(a,b), xx, w)
    end
end

struct Basis <: AbstractBasis{Float64}
    egvector::Vector{ElementGrid}
    x::Vector{Float64}
    function Basis(eg::ElementGrid...)
        #TODO add checks
        new(collect(eg), vcat(eg...))
    end
end

function Basis(a, b, ne, np)
    @assert b > a "Basis creation error - b not >a"
    step = (b-a)/ne
    eg = [ElementGrid(a+(i-1)*step, a + i*step, np)  for i in 1:ne]
    Basis(eg...)
end

function Base.show(io::IO, ::MIME"text/plain", eg::ElementGrid)
    print(io, "ElementGrid - $(length(eg.x)) points")
end

function Base.show(io::IO, ::MIME"text/plain", b::Basis)
    print(io, "Basis - $(length(b.egvector)) elements")
end

Base.size(eg::ElementGrid) = size(eg.x)
Base.size(b::Basis) = size(b.x)

Base.firstindex(b::Basis) = firstindex(b.x)
Base.firstindex(eg::ElementGrid) = firstindex(eg.x)

Base.lastindex(b::Basis) = lastindex(b.x)
Base.lastindex(eg::ElementGrid) = lastindex(eg.x)

Base.getindex(eg::ElementGrid, i::Int) = eg.x[i]
Base.getindex(b::Basis, i::Int) = b.x[i]
