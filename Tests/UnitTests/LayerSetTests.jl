# LayerSet tests

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
