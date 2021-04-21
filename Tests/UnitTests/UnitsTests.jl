@testset "Boundary Conditions" begin

    #2-vectors
    @test _2VectorInt(0,0) isa TU2VectorInt

    @test _2VectorInt(0,0) isa TU2VectorReal
    @test _2VectorFloat(0,0) isa TU2VectorReal

    @test _2VectorInt(0,0) isa TU2VectorComplex
    @test _2VectorFloat(0,0) isa TU2VectorComplex
    @test _2VectorComplex(0,0) isa TU2VectorComplex

    #3-vectors
    @test _3VectorInt(0,0,0) isa TU3VectorInt

    @test _3VectorInt(0,0,0) isa TU3VectorReal
    @test _3VectorFloat(0,0,0) isa TU3VectorReal

    @test _3VectorInt(0,0,0) isa TU3VectorComplex
    @test _3VectorFloat(0,0,0) isa TU3VectorComplex
    @test _3VectorComplex(0,0,0) isa TU3VectorComplex


    # Conversions check for right length
    @test _3VectorComplex([0,1,3]) isa _3VectorComplex
    @test_throws DimensionMismatch _3VectorComplex([0,1])
    @test_throws DimensionMismatch _3VectorComplex(_2VectorComplex(0,1))


end
