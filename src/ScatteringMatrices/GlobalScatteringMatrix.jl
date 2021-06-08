# Represents the scattering matrix of a whole stack
mutable struct GlobalScatteringMatrix{PrecisionType<:Real}

    # Size = (4*nHarmonics) x (4*nHarmonics)
    matrix::Array{Complex{PrecisionType},2}

    function GlobalScatteringMatrix{PrecisionType}(matrix::Array{Complex{PrecisionType},2}) where {PrecisionType<:Real}
        @assert isSquare(matrix)
        return new(matrix)
    end

end

function GlobalScatteringMatrix(matrix::Array{Complex{PrecisionType},2}) where {PrecisionType<:Real}
    @assert isSquare(matrix)
    return GlobalScatteringMatrix{PrecisionType}(matrix)
end


function GlobalScatteringMatrix(layerSM::LayerScatteringMatrix{PrecisionType}) where {PrecisionType<:Real}
    return GlobalScatteringMatrix{PrecisionType}(layerSM.matrix)
end

function getQuadrantSlices(sm::GlobalScatteringMatrix)
    return getQuadrantSlices(sm.matrix)
end

# Creates a empty global scattering matrix.
# Uses a regular array.  Rewrite for distributed array, etc?
# S₁₁ == S₂₂ == empty
# S₁₂ == S₂₁ == identity matrix
function initializeGlobalScatteringMatrix( PrecisionType, nHarmonics::Integer )

    Sglobal = zeros( Complex{PrecisionType}, (4*nHarmonics, 4*nHarmonics) )
    _1, _2 = getQuadrantSlices(Sglobal)
    S₁₂ = Array{Complex{PrecisionType},2}(I,(2*nHarmonics,2*nHarmonics))
    Sglobal[_1,_2] = S₁₂
    Sglobal[_2,_1] = S₁₂

    return GlobalScatteringMatrix{PrecisionType}(Sglobal)
end


# Calculate global scattering matrix for a stack of layers
# The first layer and last layers must be a SemiInfiniteLayerDefinition
function calcGlobalScatteringMatrix(simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters )

    layerStack = simulationDefinition.layerStack

    @assert first(layerStack) isa SemiInfiniteLayerDefinition
    @assert last(layerStack) isa SemiInfiniteLayerDefinition

    kVectorSet = derivedParameters.kVectorSet
    nHarmonics = numHarmonics(kVectorSet)


    # Preallocate all the data for the calculation of layer scattering matrices
    prealloc = ScatteringMatrixAllocations{simulationDefinition.PrecisionType}(numHarmonics(kVectorSet), kVectorSet)

    # Calculate values common to all layers
    prealloc.W₀ = calcW₀( numHarmonics(kVectorSet) )
    prealloc.V₀ = calcV₀( kVectorSet )
    prealloc.Sassembled.matrix = copy(calcScatteringMatrixBottom( prealloc, derivedParameters, first(layerStack), simulationDefinition.materialCollection).matrix)



    for iLayer = 2:(length(layerStack)-1)
        SnewLayer = calcScatteringMatrix(prealloc, layerStack[iLayer], simulationDefinition, derivedParameters )
        RedhefferStarProduct!( prealloc.Sassembled, SnewLayer )

    end

    Stop = calcScatteringMatrixTop( prealloc, derivedParameters, last(layerStack), simulationDefinition.materialCollection)

    RedhefferStarProduct!( prealloc.Sassembled, Stop )

    return GlobalScatteringMatrix(prealloc.Sassembled)
end



function propagateModeCoeff(Sglobal::GlobalScatteringMatrix, inputCoefficients::InputCoefficients)

    nHarmonics = numHarmonics(inputCoefficients)
    nFields = 2*nHarmonics

    sourcesCoefficients = getStackedModeCoefficients(inputCoefficients)

    scatteredCoefficients = Sglobal.matrix*sourcesCoefficients

    scatteredCoefficientsBottom = ModeCoefficientSet(scatteredCoefficients[1:nFields], BACKWARD)
    scatteredCoefficientsTop = ModeCoefficientSet(scatteredCoefficients[(nFields+1):(nFields*2)], FORWARD)

    return OutputCoefficients(scatteredCoefficientsBottom, scatteredCoefficientsTop)
end


function propagateFields( Sglobal::GlobalScatteringMatrix, inputFields::InputFields, derivedParameters::DerivedParameters )

    W₀ = calcW₀( numHarmonics(derivedParameters.kVectorSet) )
    inputCoeff = inputFields2InputCoefficients(inputFields, W₀)
    outputCoeff = propagateModeCoeff(Sglobal, inputCoeff)
    outputFields = outputCoefficients2OutputFields(outputCoeff, W₀)

    return outputFields
end
