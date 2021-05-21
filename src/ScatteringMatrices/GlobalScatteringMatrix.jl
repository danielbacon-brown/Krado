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


# TODO: Make this more memory efficient
# Calculate global scattering matrix for a stack of layers
# The first layer and last layers must be a SemiInfiniteLayerDefinition
function calcGlobalScatteringMatrix(layerStack::Vector{T}, matCol::MaterialCollection, kVectorSet::KVectorSet, gVectorSet::GvectorSet, lattice::Lattice, PrecisionType::DataType) where T<:LayerDefinition

    @assert first(layerStack) isa SemiInfiniteLayerDefinition
    @assert last(layerStack) isa SemiInfiniteLayerDefinition

    nHarmonics = numHarmonics(kVectorSet)


    # Preallocate all the data for the calculation of layer scattering matrices
    prealloc = ScatteringMatrixAllocations{PrecisionType}(numHarmonics(kVectorSet), kVectorSet)

    # Calculate values common to all layers
    # TODO: make sure these are actually used
    prealloc.W₀ = calcW₀( nHarmonics )
    prealloc.V₀ = calcV₀( kVectorSet )

    Sassembled = calcScatteringMatrixBottom( prealloc, first(layerStack), matCol, kVectorSet)


    for iLayer = 2:(length(layerStack)-1)
        SnewLayer = calcScatteringMatrix(prealloc, layerStack[iLayer], matCol, kVectorSet, gVectorSet, lattice)
        # @show SnewLayer.matrix
        RedhefferStarProduct!( Sassembled, SnewLayer )

    end
    Stop = calcScatteringMatrixTop( prealloc, last(layerStack), matCol, kVectorSet)

    RedhefferStarProduct!( Sassembled, Stop )

    return GlobalScatteringMatrix(Sassembled)
end


function calcGlobalScatteringMatrix(simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters)
    return calcGlobalScatteringMatrix(simulationDefinition.layerStack, simulationDefinition.materialCollection, derivedParameters.kVectorSet, derivedParameters.gVectorSet, simulationDefinition.lattice, simulationDefinition.PrecisionType)
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

    inputCoeff = inputFields2InputCoefficients(inputFields, derivedParameters.freeSpaceParameters.W₀)
    outputCoeff = propagateModeCoeff(Sglobal, inputCoeff)
    outputFields = outputCoefficients2OutputFields(outputCoeff, derivedParameters.freeSpaceParameters.W₀)

    return outputFields
end
