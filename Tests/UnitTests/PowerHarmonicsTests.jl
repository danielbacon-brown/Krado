@testset "Superellipse Harmonics" begin

    #Test circular 0-order:
    M, N = 0, 0
    γ = γSUPERELLIPSE["circle"]
    harmDefPow = SuperellipseHarmonicsTruncation(M,N,γ)
    mnᵢ = calcHarmonics(harmDefPow)
    # @test m_pow == [0] && n_pow == [0]
    @test mnᵢ == _2VectorInt[[0,0]]

    #Test circular +-1 order:
    M, N = 1, 1
    γ = γSUPERELLIPSE["circle"]
    harmDefPow = SuperellipseHarmonicsTruncation(M,N,γ)
    mnᵢ = calcHarmonics(harmDefPow)
    # @test m_pow == [-1, 0,0,0, 1] && n_pow == [0, -1,0,1, 0]
    @test mnᵢ == _2VectorInt[ [0,-1], [-1,0],[0,0],[1,0], [0,1], ]

end;
