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
# function calcGlobalScatteringMatrix(layerStack::LayerStack, matCol::MaterialCollection, kVectorSet::KVectorSet, gVectorSet::GvectorSet, harmonicsSet::HarmonicsSet, lattice::Lattice, PrecisionType::DataType) where T<:LayerDefinition
# function calcGlobalScatteringMatrix(layerStack::LayerStack, matCol::MaterialCollection, derivedParameters::DerivedParameters, lattice::Lattice, PrecisionType::DataType) where T<:LayerDefinition
function calcGlobalScatteringMatrix(simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters )

    layerStack = simulationDefinition.layerStack

    @assert first(layerStack) isa SemiInfiniteLayerDefinition
    @assert last(layerStack) isa SemiInfiniteLayerDefinition

    kVectorSet = derivedParameters.kVectorSet
    nHarmonics = numHarmonics(kVectorSet)


    # Preallocate all the data for the calculation of layer scattering matrices
    # prealloc = ScatteringMatrixAllocations{PrecisionType}(numHarmonics(kVectorSet), kVectorSet)
    prealloc = ScatteringMatrixAllocations{simulationDefinition.PrecisionType}(numHarmonics(kVectorSet), kVectorSet)

    # Calculate values common to all layers
    # TODO: make sure these are actually used
    # prealloc.W₀ = calcW₀( nHarmonics )
    prealloc.W₀ = calcW₀( numHarmonics(kVectorSet) )
    prealloc.V₀ = calcV₀( kVectorSet )
    # @show prealloc.W₀.matrix
    # Sassembled = calcScatteringMatrixBottom( prealloc, first(layerStack), simulationDefinition.materialCollection, kVectorSet)
    Sassembled = calcScatteringMatrixBottom( prealloc, derivedParameters, first(layerStack), simulationDefinition.materialCollection)



    for iLayer = 2:(length(layerStack)-1)
        # SnewLayer = calcScatteringMatrix(prealloc, layerStack[iLayer], matCol, derivedParameters, lattice)
        SnewLayer = calcScatteringMatrix(prealloc, layerStack[iLayer], simulationDefinition, derivedParameters )
        # @show SnewLayer.matrix
        RedhefferStarProduct!( Sassembled, SnewLayer )

    end

    # TODO: SEND THROUGH DERIVEDPARAMETERS SO IT CAN ACCESS KZNORMBOTTOM
    Stop = calcScatteringMatrixTop( prealloc, derivedParameters, last(layerStack), simulationDefinition.materialCollection)

    RedhefferStarProduct!( Sassembled, Stop )

    return GlobalScatteringMatrix(Sassembled)
end


# function calcGlobalScatteringMatrix(simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters)
#     return calcGlobalScatteringMatrix(simulationDefinition.layerStack, simulationDefinition.materialCollection, derivedParameters.kVectorSet, derivedParameters.gVectorSet, derivedParameters.harmonicsSet, simulationDefinition.lattice, simulationDefinition.PrecisionType)
# end




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
