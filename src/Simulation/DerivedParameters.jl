# Contains data that can be derived from the Simulation, but does not require calculation of scattering matrices.

mutable struct DerivedParameters

    boundaryConditions::BoundaryConditions
    harmonicsSet::HarmonicsSet
    gVectorSet::GvectorSet
    kVectorSet::KVectorSet
    freeSpaceParameters::FreeSpaceParameters

    # The z-components of k-vectors in the top and bottom layers
    # old:
    # kzBottom::Vector{ComplexF64}
    # kzTop::Vector{ComplexF64}
    # KNORM
    kzNormBottom::Vector{ComplexF64}
    kzNormTop::Vector{ComplexF64}
    kzNormGap::Vector{ComplexF64}

    # old
    # function DerivedParameters(boundaryConditions::BoundaryConditions, harmonicsSet::HarmonicsSet, gVectorSet::GvectorSet, kVectorSet::KVectorSet, freeSpaceParameters::FreeSpaceParameters, kzBottom::Vector{ComplexF64}, kzTop::Vector{ComplexF64} )
    #     return new(boundaryConditions, harmonicsSet, gVectorSet, kVectorSet, freeSpaceParameters, kzBottom, kzTop )
    # end
    # KNORM
    # function DerivedParameters(boundaryConditions::BoundaryConditions, harmonicsSet::HarmonicsSet, gVectorSet::GvectorSet, kVectorSet::KVectorSet, freeSpaceParameters::FreeSpaceParameters, kzNormBottom::Vector{ComplexF64}, kzNormTop::Vector{ComplexF64}, kzNormGap::Vector{ComplexF64} )
    function DerivedParameters(boundaryConditions::BoundaryConditions, harmonicsSet::HarmonicsSet, gVectorSet::GvectorSet, kVectorSet::KVectorSet, freeSpaceParameters::FreeSpaceParameters, kzNormBottom::Vector{<:Number}, kzNormTop::Vector{<:Number}, kzNormGap::Vector{<:Number} )
        return new(boundaryConditions, harmonicsSet, gVectorSet, kVectorSet, freeSpaceParameters, kzNormBottom, kzNormTop, kzNormGap )
    end
end

function DerivedParameters(simulationDefinition::SimulationDefinition)
    boundaryConditions = BoundaryConditions( simulationDefinition.boundaryDefinition, simulationDefinition)
    harmonicsSet = calcHarmonicsSet( simulationDefinition.harmonicsTruncation )
    gVectorSet = GvectorSet( harmonicsSet, simulationDefinition.lattice )
    kVectorSet = createKVectorSet(simulationDefinition.boundaryDefinition.wavenumber, boundaryConditions.kXY₀, simulationDefinition.boundaryDefinition.mainHarmonicOrder, gVectorSet)


    # KNORM:
    # gapLayer = UniformLayerDefinition(1,1,1)
    kzNormBottom = calckzBottom(kVectorSet, getBoundaryLayer(simulationDefinition.layerStack,BOTTOM), simulationDefinition.materialCollection, simulationDefinition.boundaryDefinition.wavenumber)
    kzNormTop = calckzTop(kVectorSet, getBoundaryLayer(simulationDefinition.layerStack,TOP), simulationDefinition.materialCollection, simulationDefinition.boundaryDefinition.wavenumber)
    kzNormGap = diag( conj(sqrt(I - kVectorSet.KxNorm.^2 - kVectorSet.KyNorm.^2)) )
    # kzNormGap = conj(sqrt(I - kVectorSet.KxNorm.^2 - kVectorSet.KyNorm.^2))

    # old:
    # kzBottom = calckz(kVectorSet, getBoundaryLayer(simulationDefinition.layerStack,BOTTOM), simulationDefinition.materialCollection, simulationDefinition.boundaryDefinition.wavenumber)
    # kzTop = calckz(kVectorSet, getBoundaryLayer(simulationDefinition.layerStack,TOP), simulationDefinition.materialCollection, simulationDefinition.boundaryDefinition.wavenumber)

    freeSpaceParameters = FreeSpaceParameters(kVectorSet)

    # KNORM:
    return DerivedParameters( boundaryConditions, harmonicsSet, gVectorSet, kVectorSet, freeSpaceParameters, kzNormBottom, kzNormTop, kzNormGap )
    # old:
    # return DerivedParameters( boundaryConditions, harmonicsSet, gVectorSet, kVectorSet, freeSpaceParameters, kzBottom, kzTop )
end



function calcInputFields( simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters )
    return calcInputFields( derivedParameters.boundaryConditions, derivedParameters.harmonicsSet, derivedParameters.kVectorSet, simulationDefinition.layerStack, simulationDefinition.materialCollection, getWavenumber(simulationDefinition) )
end



# ASSUMES USING INPUT BY ORDER BOUNDARY CONDITIONS
function getkXYZ₀(simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters)
    kzPositive = simulationDefinition.boundaryDefinition.isTop
    return getkXYZ₀(derivedParameters.boundaryConditions, simulationDefinition.layerStack, simulationDefinition.materialCollection, kzPositive)
end


function numHarmonics(derivedParameters::DerivedParameters)
    return numHarmonics(derivedParameters.harmonicsSet)
end
