@testset "Solid" begin

    matCol = MaterialCollection()
    addMaterial!(matCol, "material1", Material(ConstantPermittivity(4), ConstantPermeability(2)) )
    addMaterial!(matCol,"material2", Material(ConstantPermittivity(6)) )
    addMaterial!(matCol,"material3", Material(ConstantPermittivity(8)) )
    addMaterial!(matCol,"material4", Material(ConstantPermittivity(10)) )

    # matDict = Dict{String,AbstractMaterial}()
    # matDict["material1"] = Material(ConstantPermittivity(4), ConstantPermeability(2))
    # matDict["material2"] = Material(ConstantPermittivity(6))
    # matDict["material3"] = Material(ConstantPermittivity(8))
    # matDict["material4"] = Material(ConstantPermittivity(10))

    # matCollection = MaterialCollection(MatDict)

    λ₀ = 1.0
    wavenumber = WavenumberByλ₀(λ₀)

    circ1 = Circle([0,0], 1)
    rect1 = Rectangle([0,0],[5,1])
    poly1 = Polygon([[0,0], [2,0], [0,2]])

    circSolid = Solid(circ1, "material1")
    rectSolid = Solid(rect1, "material2", 0.5)
    polySolid = Solid(poly1, "material3", 2)

    spatialCalc = LayerPattern( [circSolid, rectSolid, polySolid], "material4" )
    @test getMaterialAtPosition( spatialCalc, [0,0]) == "material2"
    @test getMaterialAtPosition( spatialCalc, [1.5,0.25]) == "material2"
    @test getMaterialAtPosition( spatialCalc, [0,0.75]) == "material1"
    @test getMaterialAtPosition( spatialCalc, [0.9,0.9]) == "material3"
    @test getMaterialAtPosition( spatialCalc, [1.1,1.1]) == "material4"
    @test getϵμAtPosition( spatialCalc, [0,0], matCol, wavenumber)[1] == 6
    @test getϵμAtPosition( spatialCalc, [1.5,0.25], matCol, wavenumber)[1] == 6
    @test getϵμAtPosition( spatialCalc, [0,0.75], matCol, wavenumber) == (4, 2)
    @test getϵμAtPosition( spatialCalc, [0.9,0.9], matCol, wavenumber)[1] == 8
    @test getϵμAtPosition( spatialCalc, [1.1,1.1], matCol, wavenumber)[1] == 10

    # Solid using a spatially variant material
    permittivityFunc(wavenumber::Wavenumber, position) = 1 + position[X]^2 + position[Y]^2
    addMaterial!(matCol,"material5", Material( SpatialFunctionPermittivity(permittivityFunc) ) )
    # matDict["material5"] = Material( SpatialFunctionPermittivity(permittivityFunc) )
    rect2 = Rectangle([0,0],[5,5])
    rectVariantSolid = Solid(rect2, "material5", 0.5)
    spatialCalc2 = LayerPattern( rectVariantSolid, "material1" )
    @test getMaterialAtPosition( spatialCalc2, [2,2]) == "material5"
    @test getϵAtPosition( spatialCalc2, [2,2], matCol, wavenumber ) == 9
    @test getμAtPosition( spatialCalc2, [10,10], matCol, wavenumber ) == 2

end;
