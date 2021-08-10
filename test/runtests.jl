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
    @test isapprox(sum(sin.(eg) .* get_weight(eg)), 0.0; atol=1E-12)

    # Derivative accuracy
    @test all( derivative_matrix(eg) * sin.(eg) ≈ cos.(eg) )

end
