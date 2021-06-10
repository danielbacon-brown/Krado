
# abstract type LayerSetDefinition <: LayerDefinition end

# Used to create a vector of layers with blocks according to a common set of parameters
abstract type LayerSet end

mutable struct TrapezoidLayerSet <: LayerSet

    thickness::Float64
    numLayers::Int64
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
        innerMaterialName::String,
        outerMaterialName::String;
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

    layers = Vector{LayerDefinition}(undef, trapezoidLayerSet.numLayers)

    zPositions = PositionGridZ(trapezoidLayerSet.alignment, 0, trapezoidLayerSet.thickness, trapezoidLayerSet.numLayers)
    layerThickness = trapezoidLayerSet.thickness/trapezoidLayerSet.numLayers

    for iLayer in UnitRange(1,trapezoidLayerSet.numLayers)
        layerWidths = zPositions[iLayer]/trapezoidLayerSet.thickness * (trapezoidLayerSet.topLengths .- trapezoidLayerSet.bottomLengths) .+ trapezoidLayerSet.bottomLengths
        layerCenter = zPositions[iLayer]/trapezoidLayerSet.thickness * (trapezoidLayerSet.topCenter .- trapezoidLayerSet.bottomCenter) .+ trapezoidLayerSet.bottomCenter
        layers[iLayer] = PatternedLayerDefinition(trapezoidLayerSet.numDivisions, layerThickness, LayerPattern(Solid(Rectangle(layerCenter, layerWidths),trapezoidLayerSet.innerMaterialName), trapezoidLayerSet.outerMaterialName))
    end

    return layers
end



# Superellipsoid:
# (x/X)^M + (y/Y)^N + (z/Z)^O = 1
mutable struct SuperEllipsoidLayerSet <: LayerSet

    numLayers::Int64

    # Define the ellipsoid itself
    center::_3VectorFloat # center in XYZ
    ellipsoidLimits::_3VectorFloat # A,B,C in same units as center
    verticalRegionLimits::_2VectorFloat # Range of Z in which to include the dimensions
    γ::_3VectorFloat # M,N,O values

    # Define range to use


    numDivisions::_2VectorInt
    innerMaterialName::String
    outerMaterialName::String
    alignment::GridAlignment

    function SuperEllipsoidLayerSet(
        numLayers::Integer,
        # Define the ellipsoid itself
        center, # center in XYZ
        ellipsoidLimits, # A,B,C in same units as center
        verticalRegionLimits, # Range of Z in which to include the dimensions
        γ,

        numDivisions,
        innerMaterialName::String,
        outerMaterialName::String;
        alignment::GridAlignment = CENTERALIGNMENT)

        return new(
        convert(Int32,numLayers),
        # Define the ellipsoid itself
        convert(_3VectorFloat,center), # center in XYZ
        convert(_3VectorFloat,ellipsoidLimits), # A,B,C in same units as center
        convert(_2VectorFloat,verticalRegionLimits), # Range of Z in which to include the dimensions
        convert(_3VectorFloat,γ),

        convert(_2VectorInt,numDivisions),
        innerMaterialName,
        outerMaterialName,
        alignment)
    end

end

function getLayers(superEllipsoidLayerSet::SuperEllipsoidLayerSet)

    layers = Vector{LayerDefinition}(undef, superEllipsoidLayerSet.numLayers)

    zPositions = PositionGridZ(superEllipsoidLayerSet.alignment, superEllipsoidLayerSet.verticalRegionLimits[1], superEllipsoidLayerSet.verticalRegionLimits[2], superEllipsoidLayerSet.numLayers)
    layerThickness = (superEllipsoidLayerSet.verticalRegionLimits[2]-superEllipsoidLayerSet.verticalRegionLimits[1])/superEllipsoidLayerSet.numLayers

    # Superellipsoid:
    # (x/X)^M + (y/Y)^N + (z/Z)^O = 1
    # For a given z:
    # (x/X)^M + (y/Y)^N = 1 - (z/Z)^O

    for iLayer in UnitRange(1,superEllipsoidLayerSet.numLayers)
        z = zPositions[iLayer]
        zLim = 1 - (z/superEllipsoidLayerSet.ellipsoidLimits[Z])^superEllipsoidLayerSet.γ[Z]

        # Add a superellipse shape to the layer

        layerLimitX = superEllipsoidLayerSet.ellipsoidLimits[X]*zLim^(1/superEllipsoidLayerSet.γ[X])
        layerLimitY = superEllipsoidLayerSet.ellipsoidLimits[Y]*zLim^(1/superEllipsoidLayerSet.γ[Y])
        layerLimits = [layerLimitX, layerLimitY]
        center = superEllipsoidLayerSet.center[X:Y]
        layerγ = superEllipsoidLayerSet.γ[X:Y]

        superEllipse = SuperEllipse(center, layerLimits, layerγ)
        layers[iLayer] = PatternedLayerDefinition(superEllipsoidLayerSet.numDivisions, layerThickness, LayerPattern(Solid(superEllipse, superEllipsoidLayerSet.innerMaterialName), superEllipsoidLayerSet.outerMaterialName) )
    end

    return layers
end
