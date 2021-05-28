
# abstract type LayerSetDefinition <: LayerDefinition end

# Used to create a vector of layers with blocks according to a common set of parameters
abstract type LayerSet end

mutable struct TrapezoidLayerSet <: LayerSet

    thickness::Float64
    numLayers::Int32
    bottomCenter::_2VectorFloat
    topCenter::_2VectorFloat
    bottomLengths::_2VectorFloat
    topLengths::_2VectorFloat
    numDivisions::_2VectorInt
    innerMaterialName::String
    outerMaterialName::String
    alignment::GridAlignment

    function TrapezoidLayerSet(
        thickness::Real,
        numLayers::Integer,
        bottomCenter,
        topCenter,
        bottomLengths,
        topLengths,
        numDivisions,
        innerMaterialName,
        outerMaterialName;
        alignment::GridAlignment = CENTERALIGNMENT)
        return new(
            convert(Float64,thickness),
            convert(Int32,numLayers),
            _2VectorFloat(bottomCenter),
            _2VectorFloat(topCenter),
            _2VectorFloat(bottomLengths),
            _2VectorFloat(topLengths),
            _2VectorInt(numDivisions),
            innerMaterialName,
            outerMaterialName,
            alignment )
    end

end

function getLayers(trapezoidLayerSet::TrapezoidLayerSet)

    layers = Vector{LayerDefinition}(undef, numLayers)

    zPositions = PositionGridZ(trapezoidLayerSet.alignment, 0, trapezoidLayerSet.thickness, trapezoidLayerSet.numLayers)
    layerThickness = trapezoidLayerSet.thickness/numLayers

    for iLayer in UnitRange(1,numLayers)
        layerWidths = zPositions[iLayer]/trapezoidLayerSet.thickness * (trapezoidLayerSet.topLengths .- trapezoidLayerSet.bottomLengths) .+ trapezoidLayerSet.bottomLengths
        layerCenter = zPositions[iLayer]/trapezoidLayerSet.thickness * (trapezoidLayerSet.topCenter .- trapezoidLayerSet.bottomCenter) .+ trapezoidLayerSet.bottomCenter
        layers[iLayer] = PatternedLayerDefinition(numDivisions, layerThickness, LayerPattern(Solid(Rectangle(layerCenter, layerWidths),innerMaterialName), outerMaterialName))
    end

    return layers
end
