@testset "Harmonics Set" begin

    mnᵢ = [[0,0]]
    harmonicsSet = HarmonicsSet(mnᵢ)
    @test harmonicsSet.mnᵢ == [_2VectorInt(0,0)]
    @test harmonicsSet.indᵢ_mn[_2VectorInt(0,0)] == 1
    @test harmonicsSet.Δmnᵢⱼ == [_2VectorInt(0,0)]

    mnᵢ = [ [0,0],[5,3]]
    harmonicsSet = HarmonicsSet(mnᵢ)
    @test harmonicsSet.mnᵢ == [_2VectorInt(0,0), _2VectorInt(5,3)]
    @test harmonicsSet.indᵢ_mn[_2VectorInt(0,0)] == 1
    @test harmonicsSet.indᵢ_mn[_2VectorInt(5,3)] == 2
    @test harmonicsSet.Δmnᵢⱼ == [_2VectorInt(-5,-3), _2VectorInt(0,0), _2VectorInt(5,3) ]

end;
