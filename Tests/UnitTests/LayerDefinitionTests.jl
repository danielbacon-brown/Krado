@testset "LayerDefinition" begin

    # Testing calculation of coordinate grid across lattice
    divisions = [2,2]
    lattice2D = Lattice([1.0, 0], [0.5, 2])
    # grid2D = calcUniformGridPositions(lattice2D, divisions)
    grid2D = PositionGridXYleftAligned(lattice2D, divisions)
    @test grid2D.positions ≈ [ [_2VectorFloat(0.0,0.0)] [_2VectorFloat(0.25, 1.0)];
                    [_2VectorFloat(0.5,0.0)] [_2VectorFloat(0.75, 1.0)] ]
    # Testing calculation of coordinate grid across 2D lattice
    divisions1D = 4
    lattice1D = Lattice([0.4, 0.2])
    # gridX, gridY = calcUniformGridPositions(lattice1D, divisions1D)
    # grid1D = calcUniformGridPositions(lattice1D, divisions1D)
    grid1D = PositionGridXYleftAligned(lattice1D, divisions1D)
    @test grid1D.positions ≈ reshape( [ [_2VectorFloat(0.0,0.0)]; [_2VectorFloat(0.1, 0.05)];
                    [_2VectorFloat(0.2,0.1)]; [_2VectorFloat(0.3, 0.15)] ], (4,1) )
    # Testing calculation of material names across coordinate grid
    rectSolid = Solid(Rectangle([0,0],[0.5,0.5]),"material1")
    spatialPermCalc = LayerPattern([rectSolid],"material2")
    @test getMaterialAtPosition(spatialPermCalc, grid1D) == reshape(["material1"; "material1"; "material1"; "material2"], (4,1))

    # Testing calculation of ϵμ across coordinate grid
    matCol = MaterialCollection()
    # matDict = Dict{String,AbstractMaterial}()
    # matDict["material1"] = Material(ConstantPermittivity(1))
    # matDict["material2"] = Material(ConstantPermittivity(2))
    addMaterial!(matCol,"material1", Material(ConstantPermittivity(1)) )
    addMaterial!(matCol,"material2", Material(ConstantPermittivity(2)) )
    λ₀ = 1
    wavenumber = WavenumberByλ₀(λ₀)
    @test getϵμAtPosition(spatialPermCalc,grid1D,matCol,wavenumber) == reshape( [(1, 1), (1, 1), (1, 1), (2, 1)], (4,1) )
    # @test getϵμAtPosition(spatialPermCalc,grid1D,matCol,λ) == reshape( [Complex{Float64}[1, 1], [1, 1], [1, 1], [2, 1]], (4,1) )
    @test getϵAtPosition(spatialPermCalc,grid1D,matCol,wavenumber) == reshape( [Complex{Float64}(1), 1, 1, 2], (4,1) )
    @test getμAtPosition(spatialPermCalc,grid1D,matCol,wavenumber) == reshape( [Complex{Float64}(1), 1, 1, 1], (4,1) )

    # Testing the integration over G-vectors
    ϵGrid2D = getϵAtPosition(spatialPermCalc,grid2D,matCol,wavenumber)
    G2D = _2VectorFloat(0,0)
    @test calc∫xexp𝐆𝐫(ϵGrid2D,grid2D,G2D) ≈ 1.75
    G2D = _2VectorFloat(2*π,0)
    @test calc∫xexp𝐆𝐫(ϵGrid2D,grid2D,G2D) ≈ -0.25

    #Testing integration of g-vector over 1D permittivity array
    ϵGrid1D = getϵAtPosition(spatialPermCalc,grid1D,matCol,wavenumber)
    G1D = _2VectorFloat(0,0)
    @test calc∫xexp𝐆𝐫(ϵGrid1D,grid1D,G1D) ≈ 1.25

    # Integration of ϵ over g-vector set
    harmonicsSet = calcHarmonicsSet( HarmonicsTruncationByRectangle(1,1) )
    Gvectors = GvectorSet(harmonicsSet, lattice2D)
    ∫ϵexpΔ𝐆𝐫 = calc∫xexpΔ𝐆𝐫Dict(ϵGrid2D,grid2D,Gvectors)
    @test ∫ϵexpΔ𝐆𝐫[[-2,1]] ≈ -0.25

    # Calculate ⟦ϵ⟧ aka Cϵ
    harmonicsSet = calcHarmonicsSet( HarmonicsTruncationByRectangle(1,0) )
    Gvectors = GvectorSet(harmonicsSet, lattice1D)
    ∫ϵexpΔ𝐆𝐫 = calc∫xexpΔ𝐆𝐫Dict(ϵGrid1D,grid1D,Gvectors)
    Cϵᵢⱼ = assembleConvolutionMatrix( ∫ϵexpΔ𝐆𝐫, harmonicsSet )
    @test Cϵᵢⱼ ≈ Complex{Float64}[1.25 -0.25im -0.25;
        0.25im 1.25 -0.25im;
        -0.25 0.25im 1.25]

    # Calculate ⟦ϵ⟧ using more abstract functions:
    λ₀ = 1
    wavenumber = WavenumberByλ₀(λ₀)
    matCol = MaterialCollection()
    matDict = Dict{String,AbstractMaterial}()
    addMaterial!(matCol,"material1", Material(ConstantPermittivity(1)) )
    addMaterial!(matCol,"material2", Material(ConstantPermittivity(2)) )
    lattice = Lattice([0.4, 0.2])
    rectSolid = Solid(Rectangle([0,0],[0.5,0.5]),"material1")
    spatialPermCalc = LayerPattern([rectSolid],"material2")
    divisions = 4
    thickness = 2*μm
    layerDefinition = PatternedLayerDefinition(divisions, thickness, spatialPermCalc)
    harmonicsSet = calcHarmonicsSet( HarmonicsTruncationByRectangle(1,0) )
    Gvectors = GvectorSet(harmonicsSet, lattice1D)
    Cϵᵢⱼ, Cμᵢⱼ = calcConvolutionMatrices( layerDefinition, lattice, Gvectors, matCol, wavenumber )
    @test Cϵᵢⱼ ≈ Complex{Float64}[1.25 -0.25im -0.25;
        0.25im 1.25 -0.25im;
        -0.25 0.25im 1.25]

end;
