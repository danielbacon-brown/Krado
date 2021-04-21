# Need to define the abstract type before SimulationDefinition.  But the particulars need to be defined after Simulation definition
abstract type AnalysisDefinition end




# """
#     funcA(x::Int)
# I changed the docstring
# Returns double the number `x` plus `1`.
# """
# funcA(x) = 2x + 1
# funcA(x::Int) = 2x + 1




# This defines everything needed for a single simulation. e.g. measurement with a single structure, single wavelength, single AOI, etc.

# RULES: This is mutable.  Each component can be replaced without needing to change anything else.  This can make some things weird, such as BoundaryDefinition not knowing the refractive index of the semiInfinite layers.

# Should take primarily "definitions" as arguments in constructors

# Should also have "analysis" constituent that governs what kind of data is returned from this simulation when it is run
"""
    SimulationDefinition
    
SimulationDefinition(::Lattice, ::Vector{T2}, ::HarmonicsTruncation, ::BoundaryDefinition, ::MaterialCollection, ::AnalysisDefinition; PrecisionType=Float64) where T2<:LayerDefinition

SimulationDefinition is a structure that contains all of the parameters needed to define a single "experiment", e.g. a single geometry, wavelength, angle of incidence, etc.  It does not include any results or derived parameters.
"""
mutable struct SimulationDefinition
        
    lattice::Lattice
        
    layerStack::Vector{T} where T<:LayerDefinition # Stack of layers defining the structure.  First layer is reflection layer (the layer in which light is injected).  Last layer is the transmission layer.  First and last layers must be of type UniformLayerDefinition

    harmonicsTruncation::HarmonicsTruncation  # Defines how the set of harmonics to use should be defined

    boundaryDefinition::BoundaryDefinition # Parameters of boundary conditions.  Usually just a single mode (e.g. wavelength, polarization, AOI)

    materialCollection::MaterialCollection # Collection of Materials.  Materials are referred to by a string, and have a function that returns n,k.

    analysisDefinition::AnalysisDefinition # Determines the version of runSimulation() that is used to calculate the results

    PrecisionType::Type

    function SimulationDefinition(lattice::Lattice, layerStack::Vector{T2}, harmonicsTruncation::HarmonicsTruncation, boundaryDefinition::BoundaryDefinition, materialCollection::MaterialCollection, analysisDefinition::AnalysisDefinition; PrecisionType=Float64) where T2<:LayerDefinition
        return new(lattice, layerStack, harmonicsTruncation, boundaryDefinition, materialCollection, analysisDefinition, PrecisionType)
    end

end

"""
    getWavenumber( simulationDefinition )
Returns the wavenumber to be used in the simulation.
"""
function getWavenumber(simulationDefinition::SimulationDefinition)
    return simulationDefinition.boundaryDefinition.wavenumber
end

# Is this needed?
function getk₀(simulationDefinition::SimulationDefinition)
    return getk₀(getWavenumber(simulationDefinition))
end
