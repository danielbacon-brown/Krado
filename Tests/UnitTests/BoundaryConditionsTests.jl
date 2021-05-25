# BoundaryConditions tests
@testset "Boundary Definition" begin

    # Basic input mode:
    kXYZ = _3VectorFloat( 1, 0, 1 )
    k₀ = norm(kXYZ)
    θ = 1e-6
    ϕ = π/4
    A = _2VectorComplex( 1, 1im )
    n = 1
    wavenumber = WavenumberByk₀(k₀)
    mainHarmonicOrder = [0,0]
    isTop = false
    Abyϖbottom = Dict{_2VectorInt,_2VectorComplex}()
    Abyϖtop = Dict{_2VectorInt,_2VectorComplex}()
    Abyϖbottom[_2VectorInt(-1,0)] = A
    # mode1 = Mode(kXYZ, Wavenumber(k₀), A)
    # boundaryDefinition1 = InputByOrderBoundaryConditions(mode1)
    # boundaryDefinition1 = InputByOrderBoundaryDefinition(mode1, n)
    boundaryDefinition1 = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, mainHarmonicOrder, isTop, Abyϖbottom, Abyϖtop)
    @test boundaryDefinition1.θ ≈ θ
    @test boundaryDefinition1.ϕ ≈ ϕ
    # @test boundaryDefinition1.kXY₀ ≈ kXYZ[X:Y]
    @test boundaryDefinition1.Abyϖbottom[_2VectorInt(-1,0)] ≈ A

    matCol = MaterialCollection()
    addMaterial!(matCol,"Air", Material(ConstantPermittivity(1)) )

    layerBottom = SemiInfiniteLayerDefinition("Air")
    layerTop = SemiInfiniteLayerDefinition("Air")
    layerStack = [layerBottom, layerTop]

    # Test of calcModes:
    lattice1 = Lattice(2*pi)
    harmonicsDef = HarmonicsTruncationByRectangle(3,1)
    harmonicsSet = calcHarmonicsSet(harmonicsDef)
    Gvectors = GvectorSet(harmonicsSet, lattice1)
    boundaryConditions = InputByOrderBoundaryConditions( boundaryDefinition1, matCol, layerStack)
    kVectorSet = createKVectorSet(Gvectors, boundaryConditions)

    # bottomModes1, topModes1 = calcInputModes(boundaryConditions,kVectorSet,harmonicsSet)
    # @test bottomModes1[1] ≈ mode1

    # Test using multiple inputs
    amplitudeByOrderBottom = Dict{_2VectorInt,_2VectorComplex}()
    amplitudeByOrderBottom[_2VectorInt(-1,0)] = A
    amplitudeByOrderBottom[_2VectorInt(1,0)] = A
    amplitudeByOrderTop = Dict{_2VectorInt,_2VectorComplex}()
    kXYZ2 = _3VectorFloat( 0, 0, 2 )
    kXY₀2 = getXY(kXYZ2)
    k₀2 = norm(kXYZ2)
    ρ, θ, ϕ = cartesian2SphericalCoordinates(kXYZ2)
    # @test boundaryDefinition2.θ ≈ 0
    @test ϕ ≈ 0

    boundaryDefinition2 = InputByOrderBoundaryDefinition(WavenumberByk₀(k₀2), θ, ϕ, _2VectorInt(0,0), BOTTOM, amplitudeByOrderBottom, amplitudeByOrderTop)
    # @test boundaryDefinition2.kXY₀ ≈ getXY(kXYZ2)
    # @test boundaryDefinition2.θ ≈ 0
    # @test boundaryDefinition2.ϕ ≈ 2
    @test boundaryDefinition2.Abyϖbottom[_2VectorInt(1,0)] ≈ A
    @test boundaryDefinition2.Abyϖbottom[_2VectorInt(-1,0)] ≈ A

    boundaryConditions2 = InputByOrderBoundaryConditions( boundaryDefinition2, matCol, layerStack)
    kVectorSet2 = createKVectorSet(Gvectors, boundaryConditions2)
    # modesBottom2, modesTop2 = calcInputModes(boundaryConditions2, kVectorSet2, harmonicsSet)

    # @test modesBottom2[1].A ≈ A
    # @test modesBottom2[2].A ≈ A
    # NOTE: The order is not necessarily preserved when calculating modes
    # @test modesBottom2[1].kXY ≈ getkXY(kVectorSet2, getOrderIndex(harmonicsSet,_2VectorInt(-1,0)) )
    # @test modesBottom2[2].kXY ≈ getkXY(kVectorSet2, getOrderIndex(harmonicsSet,_2VectorInt(1,0)) )








end
