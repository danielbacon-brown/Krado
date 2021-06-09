@testset "Miscellaneous" begin
    # Evenly spaced values within unit cell
    @test unitLinspace(LEFTALIGNMENT, 2) ≈ [0.0, 0.5]
    @test unitLinspace(LEFTALIGNMENT, 5) ≈ [0.0, 0.2, 0.4, 0.6, 0.8]
    @test unitLinspace(CENTERALIGNMENT, 2) ≈ [0.25, 0.75]
    @test unitLinspace(CENTERALIGNMENT, 5) ≈ [0.1, 0.3, 0.5, 0.7, 0.9]
    @test unitLinspace(RIGHTALIGNMENT, 2) ≈ [0.50, 1.0]
    @test unitLinspace(RIGHTALIGNMENT, 5) ≈ [0.2, 0.4, 0.6, 0.8, 1.0]


    # Create a grid of tuples used for coordinates.  Similar to meshgrid in Matlab, but for indices.  Used with map() for iteration across arrays.
    grid1 = [ 1 2; 3 4]
    @test getGridIndices(size(grid1)) == [(1,1) (1,2); (2,1) (2,2)]

    @test bool2posNeg(true) == 1
    @test bool2posNeg(false) == -1

    v = _2VectorFloat([1,0])
    @test vectorInverse(v) ≈ [1,0]
    v = _2VectorFloat([0,2])
    @test vectorInverse(v) ≈ [0,0.5]

end;

# the S-polarization unit vector should be rotationally symetric with k.  Should not change when changing direction in z.  (SURE ABOUT THAT?)
@testset "Miscellaneous: Vector" begin

    # When k is normal vector, then take x-direction to be ŝ
    k = _3VectorFloat(0, 0, 0.5)
    @test calcŝ(k) ≈ _3VectorFloat(1,0,0)
    @test calcp̂(k) ≈ _3VectorFloat(0,1,0)
    k = _3VectorFloat(0, 0, -0.5)
    @test calcŝ(k) ≈ _3VectorFloat(1,0,0)
    @test calcp̂(k) ≈ _3VectorFloat(0,-1,0)

    # Should give approximately same results for k-vector having slight y-component
    k = _3VectorFloat(0, 1e-6, 1)
    @test calcŝ(k) ≈ _3VectorFloat(1,0,0)
    @test calcp̂(k) ≈ _3VectorFloat(0,1,-1e-6)
    k = _3VectorFloat(0, -1e-6, 1)
    @test calcŝ(k) ≈ _3VectorFloat(-1,0,0)
    @test calcp̂(k) ≈ _3VectorFloat(0,-1,-1e-6)

    # Other angles of incidence
    k = _3VectorFloat(1, 0, 2)
    @test calcŝ(k) ≈ _3VectorFloat(0,-1,0)
    @test calcp̂(k) ≈ unitize( _3VectorFloat(2,0,-1) )
    k = _3VectorFloat(-1, 0, 2)
    @test calcŝ(k) ≈ _3VectorFloat(0,1,0)
    @test calcp̂(k) ≈ unitize( _3VectorFloat(-2,0,-1) )

    # Other angles of incidence with negative kz
    k = _3VectorFloat(1, 0, -2)
    @test calcŝ(k) ≈ _3VectorFloat(0,-1,0)
    @test calcp̂(k) ≈ unitize( _3VectorFloat(-2,0,-1) )
    k = _3VectorFloat(-1, 0, -2)
    @test calcŝ(k) ≈ _3VectorFloat(0,1,0)
    @test calcp̂(k) ≈ unitize( _3VectorFloat(2,0,-1) )
    k = _3VectorFloat(-1, 0, -2)
    ŝ, p̂ = calcŝp̂(k)
    @test ŝ ≈ _3VectorFloat(0,1,0)
    @test p̂ ≈ unitize( _3VectorFloat(2,0,-1) )


    # Computing polarization vector
    k = [0,0,1]
    A = _2VectorComplex(1, 1im)
    Exyz = fieldSPtoFieldXYZ(k, A)
    @test isapprox(Exyz, [1, 1im, 0], rtol=1e-3)
    Esp = Exyz2Esp(k,Exyz)
    @test isapprox(Esp, A, rtol=1e-3)

    k2 = [0, 0.00001, 1]
    Exyz2 = fieldSPtoFieldXYZ(k2, A)
    @test isapprox(Exyz2, [1, 1im, 0], rtol=1e-3)
    Esp2 = Exyz2Esp(k2,Exyz2)
    @test isapprox(Esp2, A, rtol=1e-3)

    k3 = [0, -0.00001, 1]
    Exyz3 = fieldSPtoFieldXYZ(k3, A)
    @test isapprox(Exyz3, [-1, -1im, 0], rtol=1e-3)
    Esp3 = Exyz2Esp(k3,Exyz3)
    @test isapprox(Esp3, A, rtol=1e-3)


    # TODO: Testing complex k-vector (what are s- and p- for an evanescent wave?)
    # k = _3VectorComplex(1,1, 5*im)
    # @test false

end;


@testset "Miscellaneous: PositionGrid Conversions" begin
    lattice = Lattice([4.0,0], [0,10.0], gridAlignment = LEFTALIGNMENT)
    numDivisions = [2,2]
    # posGrid = [ [_2VectorFloat(0.0,0.0)] [_2VectorFloat(2.0,0.0)]; [_2VectorFloat(0.0,5.0)] [_2VectorFloat(2.0,5.0)] ]
    posGrid = PositionGridXY(lattice, numDivisions)
    xCoords, yCoords = linearizePositionGrid(posGrid)
    @test xCoords == [0.0, 2.0, 0.0, 2.0]
    @test yCoords == [0.0, 0, 5.0, 5.0]

    xGrid, yGrid = extractPositionGridComponents(posGrid)
    @test xGrid == [0.0 0.0; 2.0 2.0]
    @test yGrid == [0.0 5.0; 0.0 5.0]


end;


@testset "Miscellaneous: Permittivity Refractive Index Conversions" begin
    n = 1
    @test convert_n2ϵ(n) == ComplexF64(1,0)
    @test convert_ϵ2n(convert_n2ϵ(n)) == n

    n = 2
    @test convert_n2ϵ(n) == ComplexF64(4,0)
    @test convert_ϵ2n(convert_n2ϵ(n)) == n

    n = ComplexF64(2,1)
    @test convert_n2ϵ(n) == ComplexF64(3,4)
    @test convert_ϵ2n(convert_n2ϵ(n)) == n

    n = ComplexF64(2,-1)
    @test convert_n2ϵ(n) == ComplexF64(3,-4)
    @test convert_ϵ2n(convert_n2ϵ(n)) == n

    ϵ = ComplexF64(3, 4)
    μ = 1
    @test convert_ϵμ2n(ϵ, μ) == ComplexF64(2,1)
    @test convert_ϵμ2z(ϵ, μ) ≈ ComplexF64(2/5, -1/5)

    ϵ = ComplexF64(3, 4)
    μ = 4
    @test convert_ϵμ2n(ϵ, μ) == ComplexF64(4,2)
    @test convert_ϵμ2z(ϵ, μ) ≈ ComplexF64(4/5, -2/5)

    # arrays of ϵ,μ
    ϵs = [ ComplexF64(3, 4),  ComplexF64(3, 4) ]
    μs = [ 1, 4 ]
    @test convert_ϵμ2n.(ϵs, μs) ≈ [ ComplexF64(2,1), ComplexF64(4,2) ]
    @test convert_ϵμ2z.(ϵs, μs) ≈ [ ComplexF64(2/5, -1/5), ComplexF64(4/5, -2/5) ]

    # @test convert_ϵμ2z.(ϵ, μ) ≈ ComplexF64(2/5, -1/5)

end;



@testset "Miscellaneous: Spherical coordinates" begin
    kXYZ = [1, 0, 1]
    ρ, θ, ϕ = cartesian2SphericalCoordinates(kXYZ)
    @test ρ ≈ sqrt(2)
    @test θ ≈ 0
    @test ϕ ≈ π/4

    kXYZ = spherical2CartesianCoordinates(ρ, θ, ϕ)
    @test kXYZ ≈ [1, 0, 1]

    kXYZ = [0, 0, 1]
    ρ, θ, ϕ = cartesian2SphericalCoordinates(kXYZ)
    @test ρ ≈ 1
    @test θ ≈ 0
    @test ϕ ≈ 0

end;

@testset "Miscellaneous: Jones to Mueller Matrix" begin
    # horizontal polarizer
    J1 = [1 0;
        0 0]
    M1 = 0.5*[ 1 1 0 0 ;
        1 1 0 0;
        0 0 0 0;
        0 0 0 0]
    @test JonesToMuellerMatrix(J1) == M1

    # vertical polarizer
    J2 = [0 0;
        0 1]
    M2 = 0.5*[ 1 -1 0 0 ;
        -1 1 0 0;
        0 0 0 0;
        0 0 0 0]
    @test JonesToMuellerMatrix(J2) == M2

    # 45 degree polarizer
    J3 = 0.5*[1 1;
        1 1]
    M3 = 0.5*[ 1 0 1 0 ;
        0 0 0 0;
        1 0 1 0;
        0 0 0 0]
    @test JonesToMuellerMatrix(J3) == M3

    # -45 degree polarizer
    J4 = 0.5*[1 -1;
        -1 1]
    M4 = 0.5*[ 1 0 -1 0 ;
        0 0 0 0;
        -1 0 1 0;
        0 0 0 0]
    @test JonesToMuellerMatrix(J4) == M4

    # RCP - right circular polarizer
    J5 = 1/2*[1 1im;
        -1im 1]
    M5 = 0.5*[ 1 0 0 1 ;
        0 0 0 0;
        0 0 0 0;
        1 0 0 1]
    @test JonesToMuellerMatrix(J5) == M5


end;
