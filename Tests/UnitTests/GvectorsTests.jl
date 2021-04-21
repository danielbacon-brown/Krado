@testset "G-vectors" begin
    lattice = Lattice( [1,0], [0,2] )
    M, N = 1, 0
    # harmDef = HarmonicsTruncationByRectangle(M,N)
    ϖᵢ = calcHarmonicsSet(HarmonicsTruncationByRectangle(M,N))
    aGvectorsSet = GvectorSet(ϖᵢ, lattice)
    tester = [ _2VectorFloat(6.28319, 0.0), _2VectorFloat(0.0, 0.0), _2VectorFloat(-6.28319, 0.0)]
    @test isapprox(aGvectorsSet.Gᵢ, tester, atol=0.001)
end;
