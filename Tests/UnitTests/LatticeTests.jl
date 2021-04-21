@testset "Lattice" begin
    
    # Reciprocal vector calculations
    lattice1 = Lattice(_2VectorFloat(1.0,0), [0, 0.5])
    @test lattice1.G₁ ≈ [1.0,0.0]*2*pi
    @test lattice1.G₂ == [0,2]*2*pi
    lattice2 = Lattice([1.0,0], [cosd(60), sind(60)])
    @test lattice2.G₁ == [1.0, -1/(2*sind(60))]*2*pi
    @test lattice2.G₂ == [0.0, 1/(sind(60))]*2*pi

    lattice3 = Lattice( [2,0] ) # 1D lattice
    @test lattice3.G₁ ≈ [0.5,0.0]*2*pi
    @test lattice3.G₂ == [0,0]*2*pi

    lattice4 = Lattice( [0,1] )
    @test lattice4.G₁ == 2*pi*[0.0,1.0]
    @test lattice4.G₂ == [0.0,0.0]

    # Coordinate conversion
    @test convertUVtoXY(lattice1, [1,1]) ≈ [1,0.5]
    @test convertXYtoUV(lattice1, [1,0.5]) ≈ [1,1]
    @test convertUVtoXY(lattice2, [1,1]) ≈ [1+cosd(60),sind(60)]
    @test convertXYtoUV(lattice2, [1+cosd(60),sind(60)]) ≈ [1,1]

    @test convertUVtoXY(lattice3, [0.5,0]) ≈ [1,0]
    @test convertUVtoXY(lattice3, [0.5,2]) ≈ [1,0]
    @test convertXYtoUV(lattice3, [1,0]) ≈ [0.5,0]

    @test convertUVtoXY(lattice4, [2,0]) ≈ [0,2]
    @test convertUVtoXY(lattice4, [2,1]) ≈ [0,2]
    @test convertXYtoUV(lattice4, [0,2]) ≈ [2,0]

    # Coordinate offsets
    lattice5 = Lattice(_2VectorFloat(1.0,0), [0, 0.5]; originOffsetXY = [0.3,0])
    @test convertUVtoXY(lattice5, [1,1]) ≈ [1.3,0.5]
    @test convertXYtoUV(lattice5, [1.3,0.5]) ≈ [1,1]

end;
