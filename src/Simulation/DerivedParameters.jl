export DerivedParameters
# Contains data that can be derived from the Simulation, but does not require calculation of scattering matrices.

mutable struct DerivedParameters

    boundaryConditions::BoundaryConditions
    harmonicsSet::HarmonicsSet
    gVectorSet::GvectorSet
    kVectorSet::KVectorSet

    # The z-components of k-vectors in the top and bottom layers
    kzNormBottom::Vector{ComplexF64}
    kzNormTop::Vector{ComplexF64}
    kzNormGap::Vector{ComplexF64}

    function DerivedParameters(boundaryConditions::BoundaryConditions, harmonicsSet::HarmonicsSet, gVectorSet::GvectorSet, kVectorSet::KVectorSet, kzNormBottom::Vector{<:Number}, kzNormTop::Vector{<:Number}, kzNormGap::Vector{<:Number} )
        return new(boundaryConditions, harmonicsSet, gVectorSet, kVectorSet, kzNormBottom, kzNormTop, kzNormGap )
    end
end

function DerivedParameters(simulationDefinition::SimulationDefinition)
    boundaryConditions = BoundaryConditions( simulationDefinition.boundaryDefinition, simulationDefinition)
    harmonicsSet = calcHarmonicsSet( simulationDefinition.harmonicsTruncation )
    gVectorSet = GvectorSet( harmonicsSet, simulationDefinition.lattice )
    kVectorSet = createKVectorSet(simulationDefinition.boundaryDefinition.wavenumber, boundaryConditions.kXY₀, simulationDefinition.boundaryDefinition.mainHarmonicOrder, gVectorSet, harmonicsSet)


    kzNormBottom = calckzBottom(kVectorSet, getBoundaryLayer(simulationDefinition.layerStack,BOTTOM), simulationDefinition.materialCollection, simulationDefinition.boundaryDefinition.wavenumber)
    kzNormTop = calckzTop(kVectorSet, getBoundaryLayer(simulationDefinition.layerStack,TOP), simulationDefinition.materialCollection, simulationDefinition.boundaryDefinition.wavenumber)
    kzNormGap = diag( conj(sqrt(I - kVectorSet.KxNorm.^2 - kVectorSet.KyNorm.^2)) )
    return DerivedParameters( boundaryConditions, harmonicsSet, gVectorSet, kVectorSet, kzNormBottom, kzNormTop, kzNormGap )
end



function calcInputFields( simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters )
    return calcInputFields( derivedParameters.boundaryConditions, derivedParameters.harmonicsSet, derivedParameters.kVectorSet, simulationDefinition.layerStack, simulationDefinition.materialCollection, getWavenumber(simulationDefinition) )
end
function calcInputFields( simulationDefinition::SimulationDefinition )
    return calcInputFields( simulationDefinition, DerivedParameters(simulationDefinition))
end

# ASSUMES USING INPUT BY ORDER BOUNDARY CONDITIONS
function getkXYZ₀(simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters)
    kzPositive = simulationDefinition.boundaryDefinition.isTop
    return getkXYZ₀(derivedParameters.boundaryConditions, simulationDefinition.layerStack, simulationDefinition.materialCollection, kzPositive)
end


function numHarmonics(derivedParameters::DerivedParameters)
    return numHarmonics(derivedParameters.harmonicsSet)
end
