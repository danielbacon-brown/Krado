# Contains data that can be derived from the Simulation, but does not require calculation of scattering matrices.

mutable struct DerivedParameters

    boundaryConditions::BoundaryConditions
    harmonicsSet::HarmonicsSet
    gVectorSet::GvectorSet
    kVectorSet::KVectorSet
    freeSpaceParameters::FreeSpaceParameters
    
    # The z-components of k-vectors in the top and bottom layers
    kzBottom::Vector{ComplexF64}
    kzTop::Vector{ComplexF64}
    
    function DerivedParameters(boundaryConditions::BoundaryConditions, harmonicsSet::HarmonicsSet, gVectorSet::GvectorSet, kVectorSet::KVectorSet, freeSpaceParameters::FreeSpaceParameters, kzBottom::Vector{ComplexF64}, kzTop::Vector{ComplexF64} )
        return new(boundaryConditions, harmonicsSet, gVectorSet, kVectorSet, freeSpaceParameters, kzBottom, kzTop )
    end    
end

function DerivedParameters(simulationDefinition::SimulationDefinition)
    boundaryConditions = BoundaryConditions( simulationDefinition.boundaryDefinition, simulationDefinition)
    harmonicsSet = calcHarmonicsSet( simulationDefinition.harmonicsTruncation )
    gVectorSet = GvectorSet( harmonicsSet, simulationDefinition.lattice )
    kVectorSet = createKVectorSet(simulationDefinition.boundaryDefinition.wavenumber, boundaryConditions.kXY₀, simulationDefinition.boundaryDefinition.mainHarmonicOrder, gVectorSet)
    
    kzBottom = calckz(kVectorSet, getBoundaryLayer(simulationDefinition.layerStack,BOTTOM), simulationDefinition.materialCollection, simulationDefinition.boundaryDefinition.wavenumber)
    kzTop = calckz(kVectorSet, getBoundaryLayer(simulationDefinition.layerStack,TOP), simulationDefinition.materialCollection, simulationDefinition.boundaryDefinition.wavenumber)
    
    freeSpaceParameters = FreeSpaceParameters(kVectorSet)
    
    return DerivedParameters( boundaryConditions, harmonicsSet, gVectorSet, kVectorSet, freeSpaceParameters, kzBottom, kzTop, )
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
