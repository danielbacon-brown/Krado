
export runSimulation, ZeroOrderModesAnalysisDefinition, JonesMatrixAnalysisDefinition, AllModesAnalysisDefinition, RelativeReflectedTransmittedOrdersAnalysisDefinitionTransmittanceReflectanceAnalysisDefinition, PerformanceAnalysisDefinition

# Runs the simulation and returns results.  Resullts is a namedTuple.

"""
    runSimulation(simulationDefinition::SimulationDefinition)

Runs the simulation
"""
function runSimulation(simulationDefinition::SimulationDefinition)
    return runSimulation(simulationDefinition.analysisDefinition, simulationDefinition)
end


#  Calculates the modes in both input and output, on both sides
mutable struct ZeroOrderModesAnalysisDefinition <: AnalysisDefinition
    direction::Bool
    function ZeroOrderModesAnalysisDefinition( direction::Bool)
        return new(direction)
    end
end

function runSimulation(analysisDefinition::ZeroOrderModesAnalysisDefinition, simulationDefinition::SimulationDefinition)

    @assert isa(simulationDefinition.boundaryDefinition, InputByOrderBoundaryDefinition)

     # Calc common parameters
     derivedParameters = DerivedParameters(simulationDefinition)

     # Calc input fields
     inputFields = calcInputFields(simulationDefinition, derivedParameters)

     # Calc global scattering matrix
     Sglobal = calcGlobalScatteringMatrix(simulationDefinition, derivedParameters)

     # Propagate fields
     outputFields = propagateFields( Sglobal, inputFields, derivedParameters )

     # Get output fields data
     outputFieldsXYZbottom, outputFieldsXYZtop, outputFieldsSPbottom, outputFieldsSPtop = analyzeOutputFields(outputFields, derivedParameters, simulationDefinition.layerStack, simulationDefinition.materialCollection, getWavenumber(simulationDefinition) )
     # Get input fields data
     inputFieldsXYZbottom, inputFieldsXYZtop, inputFieldsSPbottom, inputFieldsSPtop = analyzeInputFields(inputFields, derivedParameters, simulationDefinition.layerStack, simulationDefinition.materialCollection, getWavenumber(simulationDefinition) )

     # Get Zero order index
     zeroOrderIndex = getOrderIndex(derivedParameters.harmonicsSet, _2VectorInt(0,0))

     if analysisDefinition.direction == FORWARD
         Tsp = outputFieldsSPtop.fields[zeroOrderIndex,:]
         Rsp = outputFieldsSPbottom.fields[zeroOrderIndex,:]
     else
         Tsp = outputFieldsSPbottom.fields[zeroOrderIndex,:]
         Rsp = outputFieldsSPtop.fields[zeroOrderIndex,:]
     end

     # Output data as a named tuple
     data = (Tsp = Tsp, Rsp = Rsp)

     return data
end

# Assumed to be zero order only
mutable struct JonesMatrixAnalysisDefinition <: AnalysisDefinition
    direction::Bool
    function JonesMatrixAnalysisDefinition( direction::Bool)
        return new(direction)
    end
end

function runSimulation(analysisDefinition::JonesMatrixAnalysisDefinition, simulationDefinition::SimulationDefinition)
# Switches the polarization of the modes

    Spol = [1, 0]
    Ppol = [0, 1]

    # assuming InputByOrderBoundaryDefinition
    # Set polarization to S:
    if analysisDefinition.direction==FORWARD
        simulationDefinition.boundaryDefinition.Abyϖbottom = Dict{_2VectorInt, _2VectorComplex}( simulationDefinition.boundaryDefinition.mainHarmonicOrder => Spol)
        simulationDefinition.boundaryDefinition.Abyϖtop = Dict{_2VectorInt, _2VectorComplex}()
    else
        simulationDefinition.boundaryDefinition.Abyϖtop = Dict{_2VectorInt, _2VectorComplex}( simulationDefinition.boundaryDefinition.mainHarmonicOrder => Spol)
        simulationDefinition.boundaryDefinition.Abyϖbottom = Dict{_2VectorInt, _2VectorComplex}()
    end

    dataS = runSimulation(ZeroOrderModesAnalysisDefinition(analysisDefinition.direction), simulationDefinition)


    # Set polarization to P:
    if analysisDefinition.direction==FORWARD
        simulationDefinition.boundaryDefinition.Abyϖbottom = Dict{_2VectorInt, _2VectorComplex}( simulationDefinition.boundaryDefinition.mainHarmonicOrder => Ppol)
        simulationDefinition.boundaryDefinition.Abyϖtop = Dict{_2VectorInt, _2VectorComplex}()
    else
        simulationDefinition.boundaryDefinition.Abyϖtop = Dict{_2VectorInt, _2VectorComplex}( simulationDefinition.boundaryDefinition.mainHarmonicOrder => Ppol)
        simulationDefinition.boundaryDefinition.Abyϖbottom = Dict{_2VectorInt, _2VectorComplex}()
    end
    dataP = runSimulation(ZeroOrderModesAnalysisDefinition(analysisDefinition.direction), simulationDefinition::SimulationDefinition)

    # Convert to Jones matrix
    # [Eₛ]  = [Mₛₛ, Mₛₚ] * [Eₛ]
    #  Eₚ      Mₚₛ, Mₚₚ     Ep

    Sinput = 1
    Pinput = 1

    # Es2 = Mss*Es1 + Msp*Ep1
    # For Spol, Ep1==0
    # Mss = Es2/Es1

    Mss = dataS.Tsp[S] / Sinput
    Msp = dataP.Tsp[S] / Pinput
    Mps = dataS.Tsp[P] / Sinput
    Mpp = dataP.Tsp[P] / Pinput

    Mtransmission = [ Mss Msp; Mps Mpp]

    # Reflection:
    Mss = dataS.Rsp[S] / Sinput
    Msp = dataP.Rsp[S] / Pinput
    Mps = dataS.Rsp[P] / Sinput
    Mpp = dataP.Rsp[P] / Pinput

    Mreflection = [ Mss Msp; Mps Mpp]

    data = ( reflection = Mreflection, transmission = Mtransmission )
    return data
end



#  Calculates the modes in both input and output, on both sides
mutable struct AllModesAnalysisDefinition <: AnalysisDefinition
    # Takes no parameters
end

function runSimulation(analysisDefinition::AllModesAnalysisDefinition, simulationDefinition::SimulationDefinition)

     # Calc common parameters
     derivedParameters = DerivedParameters(simulationDefinition)

     # Calc input fields
     inputFields = calcInputFields(simulationDefinition, derivedParameters)

     # Calc global scattering matrix
     Sglobal = calcGlobalScatteringMatrix(simulationDefinition, derivedParameters)

     # Propagate fields
     outputFields = propagateFields( Sglobal, inputFields, derivedParameters )

     # Output data as a named tuple
     data = (inputFields = inputFields, outputFields = outputFields)

     return data
 end




#  Calculates the modefields forward and backward in the top and bottom layers.  Relative to the input power.
mutable struct RelativeReflectedTransmittedOrdersAnalysisDefinition <: AnalysisDefinition
    # Takes no parameters
end

function runSimulation(analysisDefinition::RelativeReflectedTransmittedOrdersAnalysisDefinition, simulationDefinition::SimulationDefinition)

     # Calc common parameters
     derivedParameters = DerivedParameters(simulationDefinition)

     # Calc input fields
     inputFields = calcInputFields(simulationDefinition, derivedParameters)

     # Calc global scattering matrix
     Sglobal = calcGlobalScatteringMatrix(simulationDefinition, derivedParameters)

     # Propagate fields
     outputFields = propagateFields( Sglobal, inputFields, derivedParameters )

     # Get absolute flux through top and bottom layers
     inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux = calcPowerFluxes(inputFields::InputFields, outputFields::OutputFields, derivedParameters::DerivedParameters)

     # Get flux relative to input power
     inputBottomRelativeFlux, inputTopRelativeFlux, outputBottomRelativeFlux, outputTopRelativeFlux = calcRelativeFluxes(inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux)

     # Output data as a named tuple
     data = (inputBottomRelativeFlux = inputBottomRelativeFlux, inputTopRelativeFlux = inputTopRelativeFlux, outputBottomRelativeFlux = outputBottomRelativeFlux, outputTopRelativeFlux = outputTopRelativeFlux)

     return data
 end

#  Calculates the total reflectance and transmittance.
mutable struct TransmittanceReflectanceAnalysisDefinition <: AnalysisDefinition
    isForward::Bool
    function TransmittanceReflectanceAnalysisDefinition(isForward::Bool)
        return new(isForward)
    end
end

function runSimulation(analysisDefinition::TransmittanceReflectanceAnalysisDefinition, simulationDefinition::SimulationDefinition)

     # Calc common parameters
     derivedParameters = DerivedParameters(simulationDefinition)

     # Calc input fields
     inputFields = calcInputFields(simulationDefinition, derivedParameters)

     # Calc global scattering matrix
     Sglobal = calcGlobalScatteringMatrix(simulationDefinition, derivedParameters)

     # Propagate fields
     outputFields = propagateFields( Sglobal, inputFields, derivedParameters )


     # Get absolute flux through top and bottom layers
     inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux = calcPowerFluxes(inputFields::InputFields, outputFields::OutputFields, derivedParameters::DerivedParameters)

     # Get flux relative to input power
     inputBottomRelativeFlux, inputTopRelativeFlux, outputBottomRelativeFlux, outputTopRelativeFlux = calcRelativeFluxes(inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux)

     if analysisDefinition.isForward
         totalTransmittance = sum(real(outputTopRelativeFlux))
         totalReflectance = sum(real(outputBottomRelativeFlux))
     else
         totalReflectance = sum(real(outputTopRelativeFlux))
         totalTransmittance = sum(real(outputBottomRelativeFlux))
     end
     totalAbsorbance = 1 - totalReflectance - totalTransmittance


     # Output data as a named tuple
     data = (totalTransmittance = totalTransmittance, totalReflectance = totalReflectance, totalAbsorbance = totalAbsorbance)

     return data
 end



 mutable struct PerformanceAnalysisDefinition <: AnalysisDefinition

 end

 function runSimulation(analysisDefinition::PerformanceAnalysisDefinition, simulationDefinition::SimulationDefinition)

    # timerOutput = TimerOutput()

    # Calc common parameters
    derivedParameters = DerivedParameters(simulationDefinition)
    # @timeit timerOutput "derivedParameters" derivedParameters = DerivedParameters(simulationDefinition)

    # Calc input fields
    inputFields = calcInputFields(simulationDefinition, derivedParameters)

    # Calc global scattering matrix
    Sglobal = calcGlobalScatteringMatrix(simulationDefinition, derivedParameters)
    # Sglobal = calcGlobalScatteringMatrixTimed(simulationDefinition.layerStack, simulationDefinition.materialCollection, derivedParameters.kVectorSet, derivedParameters.gVectorSet, simulationDefinition.lattice)
    # Sglobal = calcGlobalScatteringMatrix(simulationDefinition.layerStack, simulationDefinition.materialCollection, derivedParameters.kVectorSet, derivedParameters.gVectorSet, simulationDefinition.lattice)



    # Propagate fields
    # @timeit timerOutput "outputFields" outputFields = propagateFields( Sglobal, inputFields, derivedParameters )
    outputFields = propagateFields( Sglobal, inputFields, derivedParameters )

    # Get absolute flux through top and bottom layers
    inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux = calcPowerFluxes(inputFields::InputFields, outputFields::OutputFields, derivedParameters::DerivedParameters)

    # Get flux relative to input power
    inputBottomRelativeFlux, inputTopRelativeFlux, outputBottomRelativeFlux, outputTopRelativeFlux = calcRelativeFluxes(inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux)

    # For verification of results
    totalTransmittance = sum(real(outputTopRelativeFlux))
    totalReflectance = sum(real(outputBottomRelativeFlux))
    totalAbsorbance = 1 - totalReflectance - totalTransmittance

    # Output data as a named tuple
    data = (totalTransmittance = totalTransmittance, totalReflectance = totalReflectance, totalAbsorbance = totalAbsorbance)

    return data
end
