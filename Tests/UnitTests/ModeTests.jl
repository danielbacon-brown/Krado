# NOT USING "MODE" ANYMORE

@testset "Mode Tests" begin

    𝐤xyz = _3VectorFloat( 1, 0, 1 )
    k₀ = Wavenumber(norm(𝐤xyz))
    𝐀 = _2VectorComplex( 1, 1im )
    mode1 = Mode(𝐤xyz, k₀,  𝐀)
    n = 1
    @test getkz(mode1, n) ≈ 1


end
