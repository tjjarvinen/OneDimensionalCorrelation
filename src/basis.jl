

abstract type AbstractElementGrid{T} <: AbstractVector{T} end
abstract type AbstractBasis{T} <: AbstractVector{T} end


struct Element1D
    a::Float64
    b::Float64
    function Element1D(a,b)
        @argcheck a < b
        new(a,b)
    end
end


struct ElementGrid{T} <: AbstractElementGrid{T}
    el::Element1D
    gl::GaussLegendre{T}
    scaling::T
    shift::T
    function ElementGrid(a, b, n)
        new{Float64}(Element1D(a,b), GaussLegendre(n), (b-a)/2, (a+b)/2)
    end
end

struct ElementGridLobatto{T} <: AbstractElementGrid{T}
    el::Element1D
    gl::LobattoLegendre{T}
    scaling::T
    shift::T
    function ElementGridLobatto(a, b, n)
        new{Float64}(Element1D(a,b), LobattoLegendre(n), (b-a)/2, (a+b)/2)
    end
end


struct Basis{T} <: AbstractBasis{T}
    egvector::Vector{ElementGrid{T}}
    x::Vector{T}
    function Basis(eg::ElementGrid{T}...) where T
        if length(eg) > 1
            for i in 2:length(eg)
                @assert eg[i-1] <= eg[i]
            end
        end
        new{T}(collect(eg), vcat(eg...))
    end
end


struct BasisLobatto{T} <: AbstractBasis{T}
    egvector::Vector{ElementGridLobatto{T}}
    x::Vector{T}
    D::Matrix{T}
    w::Vector{T}
    function BasisLobatto(eg::ElementGridLobatto{T}...) where T
        if length(eg) > 1
            x = Vector(eg[1])
            for i in 2:length(eg)
                @assert eg[i-1][end] ≈ eg[i][begin]
                append!(x, eg[i][begin+1:end])
            end

            # Build derivative matrix and integration weights
            dv = [ derivative_matrix(x) for x in eg ]
            cr = [size(x,1) for x in dv] # dv is square matrix
            # Size of derivative matrix. Lobatto basis overlaps with adjacent elements
            # and we need to remove this overlap, by reducing number of points.
            nd = sum(cr) - ( length(cr) - 1 )
            D = zeros(T ,nd, nd)  # Derivate matrix
            w = zeros(T, nd)      # Integration weights
            D[1:cr[1], 1:cr[1]] = dv[1]
            w[1:cr[1]] = get_weight(eg[1])
            # Index range for final derivative matrix build from element ones
            tmp = ones(Int, length(cr))
            tmp[1] = 0
            cum_index = cumsum( cr .- tmp )
            for i in 2:length(cum_index)
                ir = cum_index[i-1]:cum_index[i]
                D[ir, ir] += dv[i]
                D[cum_index[i-1],:] .*= 0.5 # remove double count from overlap point
                w[ir] += get_weight(eg[i])
            end
        else
            x = Vector(eg[1])
            D = eg[1].gl.D
            w = get_weight(eg[1])
        end
        new{T}(collect(eg), x, D, w)
    end
end


function Basis(a, b, ne, np)
    @argcheck a < b 
    step = (b-a)/ne
    eg = [ElementGrid( a+(i-1)*step, a + i*step, np )  for i in 1:ne]
    return Basis(eg...)
end

function BasisLobatto(a, b, ne, np)
    @argcheck a < b 
    step = (b-a)/ne
    eg = [ElementGridLobatto( a+(i-1)*step, a + i*step, np )  for i in 1:ne]
    return BasisLobatto(eg...)
end


function Base.show(io::IO, ::MIME"text/plain", eg::ElementGrid)
    print(io, "ElementGrid - $(length(eg)) points")
end

function Base.show(io::IO, ::MIME"text/plain", b::Basis)
    print(io, "Basis - $(length(b.egvector)) elements and $(length(b)) points")
end


Base.size(eg::AbstractElementGrid) = size(eg.gl.nodes)
Base.size(b::AbstractBasis) = size(b.x)

Base.firstindex(b::AbstractBasis) = firstindex(b.x)
Base.firstindex(eg::AbstractElementGrid) = firstindex(eg.gl.nodes)

Base.lastindex(b::AbstractBasis) = lastindex(b.x)
Base.lastindex(eg::AbstractElementGrid) = lastindex(eg.gl.nodes)

Base.getindex(eg::AbstractElementGrid, i::Int) = eg.gl.nodes[i]*eg.scaling + eg.shift
Base.getindex(b::AbstractBasis, i::Int) = b.x[i]


for OP in [:(<), :(<=)]
    @eval begin
        function Base.$OP(a::Element1D, b::Element1D)
            return $OP(a.b, b.a)
        end
        function Base.$OP(a::AbstractElementGrid, b::AbstractElementGrid)
            return $OP(a.el, b.el)
        end
    end
end

for OP in [:(>), :(>=)]
    @eval begin
        function Base.$OP(a::Element1D, b::Element1D)
            return $OP(a.a, b.b)
        end
        function Base.$OP(a::AbstractElementGrid, b::AbstractElementGrid)
            return $OP(a.el, b.el)
        end
    end
end


function get_weight(eg::AbstractElementGrid, i::Int)
    return eg.gl.weights[i] * eg.scaling
end

function get_weight(eg::AbstractElementGrid)
    return eg.gl.weights .* eg.scaling
end

function PolynomialBases.derivative_matrix(eg::AbstractElementGrid)
    return eg.gl.D ./ eg.scaling
end

function metric_tensor(b::AbstractBasis)
    return Diagonal(get_weight(b))
end

function PolynomialBases.derivative_matrix(b::Basis)
    dv = [ derivative_matrix(x) for x in b.egvector  ]
    # get block sizes for block banded matrix
    cr = [size(x)[1] for x in dv] # dv is square matrix
    out = BlockBandedMatrix{eltype(b)}(undef, cr, cr, (0,0))
    for i in 1:length(cr)
        out[Block(i,i)] = dv[i]
    end
    return Matrix(out) 
end


function PolynomialBases.derivative_matrix(b::BasisLobatto)
    return b.D
end


"""
    get_weight(b::AbstractBasis)

Integration weights for the basis
"""
function get_weight(b::Basis)
    tmp = get_weight.(b.egvector)
    return vcat(tmp...)
end


function get_weight(b::BasisLobatto)
    return b.w
end

"""
    get_identity(b::AbstractBasis)

Identity matrix for the given basis
"""
function get_identity(b::AbstractBasis)
    return diagm(ones(length(b)))
end


"""
    get_length(el::Element1D)
    get_length(eg::ElementGrid)
    get_length(b::Basis)

Gives physical length the input represents
"""
function get_length(el::Element1D)
    return el.b - el.a
end

function get_length(eg::AbstractElementGrid)
    return get_length(eg.el)
end

function get_length(b::AbstractBasis)
    return sum( get_length, b.egvector) 
end