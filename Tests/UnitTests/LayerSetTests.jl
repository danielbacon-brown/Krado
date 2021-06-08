# LayerSet tests

@testset "TrapezoidLayerSet" begin

thickness = 0.9
numLayers = 9
bottomCenter = [0.5, 0.5]
topCenter = [-0.5, -0.5]
bottomLengths = [ 1.0, 2.0]
topLengths = [1.5, 2.5]
numDivisions = [100,100]
innerMaterialName = "Si"
outerMaterialName = "Air"
trapezoidalLayerSet = TrapezoidLayerSet(
    thickness,
    numLayers,
    bottomCenter,
    topCenter,
    bottomLengths,
    topLengths,
    numDivisions,
    innerMaterialName,
    outerMaterialName;
    alignment=CENTERALIGNMENT)

layers = getLayers(trapezoidalLayerSet)
middleLayer = layers[5]
@test middleLayer.thickness ≈ .1
@test middleLayer.layerPattern.solids[1].shape.center == [0,0]
@test middleLayer.layerPattern.solids[1].shape.lengths ≈ [1.25,2.25]
@test getMaterialAtPosition( middleLayer.layerPattern, [0.6, 1.1]) == "Si"
@test getMaterialAtPosition( middleLayer.layerPattern, [0.7, 1.1]) == "Air"
end;



# Superellipse layer test sets
@testset "SuperEllipsoidLayerSet" begin
    numLayers=5
    center=[0,0,0] # center in XYZ
    ellipsoidLimits=[1,2,4] # A,B,C in same units as center
    verticalRegionLimits=[-1,1] # Range of Z in which to include the dimensions
    γ=[1,2,3]
    numDivisions=[100,100]
    innerMaterialName="Si"
    outerMaterialName="Air"
    superEllipsoid = SuperEllipsoidLayerSet(
        numLayers::Integer,
        # Define the ellipsoid itself
        center, # center in XYZ
        ellipsoidLimits, # A,B,C in same units as center
        verticalRegionLimits, # Range of Z in which to include the dimensions
        γ,
        numDivisions,
        innerMaterialName,
        outerMaterialName)
    layers = getLayers(superEllipsoid)
    middleLayer = layers[5]
    @test middleLayer.thickness ≈ .4
    @test getMaterialAtPosition( middleLayer.layerPattern, [0, 0]) == "Si"
    @test getMaterialAtPosition( middleLayer.layerPattern, [1.5, 0]) == "Air"
    @test getMaterialAtPosition( middleLayer.layerPattern, [0, 1.5]) == "Si"
    @test getMaterialAtPosition( middleLayer.layerPattern, [0.5, 1]) == "Si"

end;
