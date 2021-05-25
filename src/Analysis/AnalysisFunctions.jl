## Helper functions analyzing data

function calcPowerFluxes(inputFields::InputFields, outputFields::OutputFields, derivedParameters::DerivedParameters)

    bottomFieldsOutput = convertFieldSetStackToXYZ(outputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzNormBottom)
    topFieldsOutput = convertFieldSetStackToXYZ(outputFields.top, derivedParameters.kVectorSet, derivedParameters.kzNormTop)
    bottomFieldsInput = convertFieldSetStackToXYZ(inputFields.bottom, derivedParameters.kVectorSet, -derivedParameters.kzNormBottom)
    topFieldsInput = convertFieldSetStackToXYZ(inputFields.top, derivedParameters.kVectorSet, -derivedParameters.kzNormTop)

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
# TODO: unit test this: and combine into single function
function analyzeOutputFields(outputFields::OutputFields, derivedParameters::DerivedParameters, layerStack::Vector{<:LayerDefinition}, matCol::MaterialCollection, wavenumber::Wavenumber )

    outputFieldSetXYZbottom = convertFieldSetStackToXYZ(outputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzNormBottom)
    outputFieldSetXYZtop = convertFieldSetStackToXYZ(outputFields.top, derivedParameters.kVectorSet, derivedParameters.kzNormTop)

    outputFieldSetSPbottom = convertFieldSetXYZtoSP( outputFieldSetXYZbottom, derivedParameters.kVectorSet, getMaterial(matCol, first(layerStack).backgroundMaterialName), wavenumber )
    outputFieldSetSPtop = convertFieldSetXYZtoSP( outputFieldSetXYZtop, derivedParameters.kVectorSet, getMaterial(matCol, last(layerStack).backgroundMaterialName), wavenumber )

    return outputFieldSetXYZbottom, outputFieldSetXYZtop, outputFieldSetSPbottom, outputFieldSetSPtop
end

# clone of output fields (is it exact?)
function analyzeInputFields( inputFields::InputFields, derivedParameters::DerivedParameters, layerStack::Vector{<:LayerDefinition}, matCol::MaterialCollection, wavenumber::Wavenumber )

    # Convert xyz fields to stacked fields
    inputFieldSetXYZbottom = convertFieldSetStackToXYZ(inputFields.bottom, derivedParameters.kVectorSet, -derivedParameters.kzNormBottom)
    inputFieldSetXYZtop = convertFieldSetStackToXYZ(inputFields.top, derivedParameters.kVectorSet, -derivedParameters.kzNormTop)

    inputFieldSetSPbottom = convertFieldSetXYZtoSP( inputFieldSetXYZbottom, derivedParameters.kVectorSet, getMaterial(matCol, first(layerStack).backgroundMaterialName), wavenumber )
    inputFieldSetSPtop = convertFieldSetXYZtoSP( inputFieldSetXYZtop, derivedParameters.kVectorSet, getMaterial(matCol, last(layerStack).backgroundMaterialName), wavenumber )

    return inputFieldSetXYZbottom, inputFieldSetXYZtop, inputFieldSetSPbottom, inputFieldSetSPtop
end
