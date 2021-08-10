

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


struct ElementGrid{T} <: AbstractElementGrid{T}
    el::Element1D
    gl::GaussLegendre{T}
    shift::T
    scaling::T
    function ElementGrid(a, b, n)
        new{Float64}(Element1D(a,b), GaussLegendre(n), (b-a)/2, (a+b)/2)
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
    print(io, "ElementGrid - $(length(eg)) points")
end

function Base.show(io::IO, ::MIME"text/plain", b::Basis)
    print(io, "Basis - $(length(b)) elements")
end


Base.size(eg::ElementGrid) = size(eg.gl.nodes)
Base.size(b::Basis) = size(b.x)

Base.firstindex(b::Basis) = firstindex(b.x)
Base.firstindex(eg::ElementGrid) = firstindex(eg.gl.nodes)

Base.lastindex(b::Basis) = lastindex(b.x)
Base.lastindex(eg::ElementGrid) = lastindex(eg.gl.nodes)

Base.getindex(eg::ElementGrid, i::Int) = eg.gl.nodes[i]*eg.scaling + eg.shift
Base.getindex(b::Basis, i::Int) = b.x[i]


for OP in [:(<), :(<=)]
    @eval begin
        function Base.$OP(a::Element1D, b::Element1D)
            return $OP(a.b, b.a)
        end
        function Base.$OP(a::ElementGrid, b::ElementGrid)
            return $OP(a.el, b.el)
        end
    end
end

for OP in [:(>), :(>=)]
    @eval begin
        function Base.$OP(a::Element1D, b::Element1D)
            return $OP(a.a, b.b)
        end
        function Base.$OP(a::ElementGrid, b::ElementGrid)
            return $OP(a.el, b.el)
        end
    end
end


function get_weight(eg::ElementGrid, i::Int)
    return eg.gl.weights[i] * eg.scaling
end

function get_weight(eg::ElementGrid)
    return eg.gl.weights .* eg.scaling
end

function derivative_matrix(eg::ElementGrid)
    return eg.gl.D ./ eg.scaling
end