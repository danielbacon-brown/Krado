@testset "Rectangular Harmonics" begin

    #Test only 0-order:
    M, N = 0, 0
    harmDef = HarmonicsTruncationByRectangle(M,N)
    mnᵢ = calcHarmonics(harmDef)
    @test mnᵢ == [_2VectorInt(0,0)]

    # Test +-5 M orders:
    M, N = 5, 0
    harmDef = HarmonicsTruncationByRectangle(M,N)

    mnᵢ = calcHarmonics(harmDef)
    @test mnᵢ == _2VectorInt[[5, 0], [4, 0], [3, 0], [2, 0], [1, 0], [0, 0], [-1, 0], [-2, 0], [-3, 0], [-4, 0], [-5, 0]]

    # Test +-5 N orders:
    M, N = 0, 5
    harmDef = HarmonicsTruncationByRectangle(M,N)
    mnᵢ = calcHarmonics(harmDef)
    @test mnᵢ == _2VectorInt[[0,5], [0,4], [0,3], [0,2], [0,1], [0,0], [0,-1], [0,-2], [0,-3], [0,-4], [0,-5]]

    # Test +-1 M,N  orders
    M, N = 1, 1
    harmDef = HarmonicsTruncationByRectangle(M,N)
    mnᵢ = calcHarmonics(harmDef)
    @test mnᵢ == _2VectorInt[[1,1],[0,1],[-1,1], [1,0],[0,0],[-1,0], [1,-1],[0,-1],[-1,-1]]

end;
