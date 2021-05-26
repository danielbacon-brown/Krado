# Sweep functions

abstract type SimulationSweep end



# A sweep of simulations over a series of wavelengths
mutable struct WavenumberSweep <: SimulationSweep 

    # Iterable of wavenumbers
    wavenumbers
    
    # The simulation definition that will be modified with the given wavenumber.
    simulationBase::SimulationDefinition
    
end


function runSweep(sweep::WavenumberSweep)
    
    sweepResults = Vector{}
    # for wavenumber in sweep.wavenumbers
    #     sweep.simulationDefinition.wavenumber = wavenumber
    #     simResults = runSimulation(simulationDefinition)
    # end
    
    function doIndividualSimulation(simulationDefinition, wavenumber)
        simulationDefinition.wavenumber = wavenumber
        return runSimulation(simulationDefinition)
    end
    
    sweepResults = ( doIndividualSimulation(sweep.simulationDefinition) )
    
end




# Defines a sweep using a baseline simulationDefinition and an iterable of parameters.  Mutates the simulation definition
mutable struct MutationSimulationSweep <: SimulationSweep

    # Parameters is an iterable.  May be nested (e.g. a named tuple inside a vector).
    parameters
    
    # Baseline Simulation Sweep
    baselineSimulation::SimulationDefinition
    
    # The function that 
    mutationFunction::Function
    
end

function runSweep(sweep::MutationSimulationSweep)
    
    numElements = length(sweep.parameters)
    sweepData = Vector(UndefInitializer, numElements)
    
    for i = 1:numElements
        baselineSimulation
    end
    
end
