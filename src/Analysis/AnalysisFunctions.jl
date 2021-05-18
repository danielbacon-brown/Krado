## Helper functions analyzing data

function calcPowerFluxes(inputFields::InputFields, outputFields::OutputFields, derivedParameters::DerivedParameters)

    # Convert to 3-vector data
    #old:
    # bottomFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzBottom)
    # topFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.top, derivedParameters.kVectorSet, derivedParameters.kzTop)
    # bottomFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzBottom)
    # topFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.top, derivedParameters.kVectorSet, derivedParameters.kzTop)

    # TODO: ADJUST FOR NORM:
    # KNORM:
    bottomFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzNormBottom)
    topFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.top, derivedParameters.kVectorSet, derivedParameters.kzNormTop)
    bottomFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzNormBottom)
    topFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.top, derivedParameters.kVectorSet, derivedParameters.kzNormTop)

    # old
    # Get power flux of each
    # inputBottomPowerFlux = calcPowerFlux(bottomFieldsInput, derivedParameters.kzBottom)
    # outputBottomPowerFlux = calcPowerFlux(bottomFieldsOutput, derivedParameters.kzBottom)
    # inputTopPowerFlux = calcPowerFlux(topFieldsInput, derivedParameters.kzTop)
    # outputTopPowerFlux = calcPowerFlux(topFieldsOutput, derivedParameters.kzTop)
    # KNORM
    # Get power flux of each
    inputBottomPowerFlux = calcPowerFlux(bottomFieldsInput, derivedParameters.kzNormBottom)
    outputBottomPowerFlux = calcPowerFlux(bottomFieldsOutput, derivedParameters.kzNormBottom)
    inputTopPowerFlux = calcPowerFlux(topFieldsInput, derivedParameters.kzNormTop)
    outputTopPowerFlux = calcPowerFlux(topFieldsOutput, derivedParameters.kzNormTop)

    totalInputPowerFlux = sum(real(inputBottomPowerFlux)) + sum(real(inputTopPowerFlux))

    return inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux
end

function calcRelativeFluxes(inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux)

    totalInputPowerFlux = abs.(sum(real(inputBottomPowerFlux))) + abs.(sum(real(inputTopPowerFlux)))

    inputBottomRelativeFlux = abs.(real(inputBottomPowerFlux)) / totalInputPowerFlux
    inputTopRelativeFlux = abs.(real(inputTopPowerFlux)) / totalInputPowerFlux
    outputBottomRelativeFlux = abs.(real(outputBottomPowerFlux)) / totalInputPowerFlux
    outputTopRelativeFlux = abs.(real(outputTopPowerFlux)) / totalInputPowerFlux
    return inputBottomRelativeFlux, inputTopRelativeFlux, outputBottomRelativeFlux, outputTopRelativeFlux
end



# Returns the given outputFields in terms of XYZ and SP.
# TODO: unit test this:
function analyzeOutputFields(outputFields::OutputFields, derivedParameters::DerivedParameters, layerStack::Vector{<:LayerDefinition}, matCol::MaterialCollection, wavenumber::Wavenumber )

    # Convert xyz fields to stacked fields
    #old:
    # outputFieldSetXYZbottom = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzBottom)
    # outputFieldSetXYZtop = convertFieldSetStackToFieldSetXYZ(outputFields.top, derivedParameters.kVectorSet, derivedParameters.kzTop)
    # KNORM
    outputFieldSetXYZbottom = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzNormBottom)
    outputFieldSetXYZtop = convertFieldSetStackToFieldSetXYZ(outputFields.top, derivedParameters.kVectorSet, derivedParameters.kzNormTop)

    # outputFieldSetSPbottom, outputFieldSetSPtop = convertXYZoutputFieldsToSP(outputFieldSetXYZbottom, outputFieldSetXYZtop, derivedParameters.kVectorSet, layerStack, matCol, wavenumber)
    outputFieldSetSPbottom = convertFieldSetXYZtoSP( outputFieldSetXYZbottom, derivedParameters.kVectorSet, getMaterial(matCol, first(layerStack).backgroundMaterialName), wavenumber )
    outputFieldSetSPtop = convertFieldSetXYZtoSP( outputFieldSetXYZtop, derivedParameters.kVectorSet, getMaterial(matCol, last(layerStack).backgroundMaterialName), wavenumber )

    return outputFieldSetXYZbottom, outputFieldSetXYZtop, outputFieldSetSPbottom, outputFieldSetSPtop
end

# TODO: unit test this:
# clone of output fields (is it exact?)
function analyzeInputFields( inputFields::InputFields, derivedParameters::DerivedParameters, layerStack::Vector{<:LayerDefinition}, matCol::MaterialCollection, wavenumber::Wavenumber )

    # Convert xyz fields to stacked fields
    #old:
    # outputFieldSetXYZbottom = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzBottom)
    # outputFieldSetXYZtop = convertFieldSetStackToFieldSetXYZ(outputFields.top, derivedParameters.kVectorSet, derivedParameters.kzTop)
    # KNORM
    inputFieldSetXYZbottom = convertFieldSetStackToFieldSetXYZ(inputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzNormBottom)
    inputFieldSetXYZtop = convertFieldSetStackToFieldSetXYZ(inputFields.top, derivedParameters.kVectorSet, derivedParameters.kzNormTop)
    # @show inputFieldSetXYZbottom.fields[P]
    # outputFieldSetSPbottom, outputFieldSetSPtop = convertXYZoutputFieldsToSP(outputFieldSetXYZbottom, outputFieldSetXYZtop, derivedParameters.kVectorSet, layerStack, matCol, wavenumber)
    inputFieldSetSPbottom = convertFieldSetXYZtoSP( inputFieldSetXYZbottom, derivedParameters.kVectorSet, getMaterial(matCol, first(layerStack).backgroundMaterialName), wavenumber )
    inputFieldSetSPtop = convertFieldSetXYZtoSP( inputFieldSetXYZtop, derivedParameters.kVectorSet, getMaterial(matCol, last(layerStack).backgroundMaterialName), wavenumber )
    # @show inputFieldSetSPbottom.fields[P]

    return inputFieldSetXYZbottom, inputFieldSetXYZtop, inputFieldSetSPbottom, inputFieldSetSPtop
end
