## Helper functions analyzing data

function calcPowerFluxes(inputFields::InputFields, outputFields::OutputFields, derivedParameters::DerivedParameters)

    # Convert to 3-vector data
    bottomFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzBottom)
    topFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.top, derivedParameters.kVectorSet, derivedParameters.kzTop)
    bottomFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzBottom)
    topFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.top, derivedParameters.kVectorSet, derivedParameters.kzTop)

    # Get power flux of each
    inputBottomPowerFlux = calcPowerFlux(bottomFieldsInput, derivedParameters.kzBottom)
    outputBottomPowerFlux = calcPowerFlux(bottomFieldsOutput, derivedParameters.kzBottom)
    inputTopPowerFlux = calcPowerFlux(topFieldsInput, derivedParameters.kzTop)
    outputTopPowerFlux = calcPowerFlux(topFieldsOutput, derivedParameters.kzTop)

    totalInputPowerFlux = sum(real(inputBottomPowerFlux)) + sum(real(inputTopPowerFlux))

    return inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux 
end

function calcRelativeFluxes(inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux)

    totalInputPowerFlux = sum(real(inputBottomPowerFlux)) + sum(real(inputTopPowerFlux))

    inputBottomRelativeFlux = real(inputBottomPowerFlux) / totalInputPowerFlux
    inputTopRelativeFlux = real(inputTopPowerFlux) / totalInputPowerFlux
    outputBottomRelativeFlux = real(outputBottomPowerFlux) / totalInputPowerFlux
    outputTopRelativeFlux = real(outputTopPowerFlux) / totalInputPowerFlux
    return inputBottomRelativeFlux, inputTopRelativeFlux, outputBottomRelativeFlux, outputTopRelativeFlux
end



# Returns the given outputFields in terms of XYZ and SP.
function analyzeOutputFields(outputFields::OutputFields, derivedParameters::DerivedParameters, layerStack::Vector{<:LayerDefinition}, matCol::MaterialCollection, wavenumber::Wavenumber )

    # Convert xyz fields to stacked fields
    outputFieldSetXYZbottom = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzBottom)
    outputFieldSetXYZtop = convertFieldSetStackToFieldSetXYZ(outputFields.top, derivedParameters.kVectorSet, derivedParameters.kzTop)

    outputFieldSetSPbottom, outputFieldSetSPtop = convertXYZoutputFieldsToSP(outputFieldSetXYZbottom, outputFieldSetXYZtop, derivedParameters.kVectorSet, layerStack, matCol, wavenumber)


    return outputFieldSetXYZbottom, outputFieldSetXYZtop, outputFieldSetSPbottom, outputFieldSetSPtop
end
