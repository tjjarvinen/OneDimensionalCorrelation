using OneDimensionalCorrelation
using Test

@testset "Elements" begin
    e1 = Element1D(0, 1.2)
    e2 = Element1D(1.2, 2.4)
    
    @test (e1 < e2) == false
    @test (e2 > e1) == false
    @test e1 <= e2
    @test e2 >= e1

    p = 64
    eg1 = ElementGrid(0, 1.2, p)
    eg2 = ElementGrid(1.2, 2.4, p)

    @test (eg1 < eg2) == false
    @test (eg2 > eg1) == false
    @test eg1 <= eg2
    @test eg2 >= eg1

    eg = ElementGrid(0, 2π, p)

    # Integration accuracy
    @test isapprox( sum(sin.(eg) .* get_weight(eg)), 0.0; atol=1E-12 )

    # Derivative accuracy
    @test all( derivative_matrix(eg) * sin.(eg) ≈ cos.(eg) )

end

@testset "Basis Gauss-Lagrange" begin
    b = Basis(-π, π, 2, 64)

    # Symmetric range 
    @test b[begin] == -b[end]

    w = get_weight(b)
    ∇ = derivative_matrix(b)

    # Integral
    @test isapprox( w' *  sin.(b), 0.0; atol=1E-12 )

    # Derivative
    @test ∇ * sin.(b) ≈ cos.(b) 

end

@testset "Basis Gauss-Lobatto" begin
    b = BasisLobatto(-2, 2, 3, 32)

    # Symmetric range 
    @test b[begin] == -b[end]

    w = get_weight(b)
    ∇ = derivative_matrix(b)

    # Integral
    @test isapprox( w' *  sin.(b), 0.0; atol=1E-12 )

    # Derivative
    @test  ∇ * sin.(b) ≈ cos.(b) 

end

@testset "Integrals" begin
    b = Basis(-5, 5, 2, 64)

    ψ = exp.(-b.^2)
    I = get_identity(b)

    # Norm of Gaussian
    @test bracket(b, ψ, I, ψ) ≈ sqrt(π/2)
    @test bracket(b, ψ, ψ) ≈ sqrt(π/2)

    # Particle in box states
    ϕ1 = particle_in_box(b, 1)
    ϕ4 = particle_in_box(b, 4)
    @test bracket(b, ϕ1, ϕ1) ≈ 1
    @test bracket(b, ϕ4, ϕ4) ≈ 1
    @test bracket(b, ϕ1, ϕ4) + 1 ≈ 1

    # Particle in box Energy
    E(b::Basis, n::Int) = 0.5 * n^2 * π^2 / (get_length(b))^2

    d = derivative_matrix(b)
    e_kin = -0.5 * d^2

    # Kinetic energy accuracy
    abs( bracket(b, ϕ1, e_kin, ϕ1) - E(b, 1) ) < 1E-9
    abs( bracket(b, ϕ4, e_kin, ϕ4) - E(b, 4) ) < 1E-9
end

@testset "Hartree Fock" begin
    b = BasisLobatto(-5, 5, 2, 64)
    orbitals = initial_orbitals(b)
    C = coulomb_matrix(b, orbitals)
    K = exchange_matrix(b, orbitals)    
    F = fock_matrix(b)

    # Matrices are Hermitian
    @test C[2,1] ≈ C[2,1]
    @test K[2,1] ≈ K[2,1]
    @test F[2,1] ≈ F[2,1]  
end