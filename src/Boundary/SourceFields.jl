# Container for data describing the fields entering the device.
mutable struct InputFields

    # The mode fields in the bottom and top layers
    # bottom must be travelling forward
    # top must be travelling backward
    bottom::FieldSetStack
    top::FieldSetStack

    function InputFields(bottom::FieldSetStack, top::FieldSetStack )

        # Make sure the modesets are in the right direction
        @assert bottom.isForward == true
        @assert top.isForward == false

        return new(bottom, top)
    end
end

# Container for data describing the fields leaving the device.
mutable struct OutputFields

    # The mode fields in the bottom and top layers
    # bottom must be travelling backward
    # top must be travelling forward
    bottom::FieldSetStack
    top::FieldSetStack

    function OutputFields(bottom::FieldSetStack, top::FieldSetStack )

        # Make sure the modesets are in the right direction
        @assert bottom.isForward == false
        @assert top.isForward == true

        return new(bottom, top)
    end
end


# Container for data describing the coefficients entering the device
mutable struct InputCoefficients

    # The mode fields in the bottom and top layers
    # bottom must be travelling forward
    # top must be travelling backward
    bottom::ModeCoefficientSet
    top::ModeCoefficientSet

    function InputCoefficients(bottom::ModeCoefficientSet, top::ModeCoefficientSet )

        # Make sure the modesets are in the right direction
        @assert bottom.isForward == true
        @assert top.isForward == false

        @assert length(bottom.modeCoefficients) == length(top.modeCoefficients)

        return new(bottom, top)
    end
end

# Container for data describing the coefficients leaving the device
mutable struct OutputCoefficients

    # The mode fields in the bottom and top layers
    # bottom must be travelling backward
    # top must be travelling forward
    bottom::ModeCoefficientSet
    top::ModeCoefficientSet

    function OutputCoefficients(bottom::ModeCoefficientSet, top::ModeCoefficientSet )

        # Make sure the modesets are in the right direction
        @assert bottom.isForward == false
        @assert top.isForward == true

        return new(bottom, top)
    end
end


function getStackedModeCoefficients( inputCoefficients::InputCoefficients)::Vector{ComplexF64}
    return vcat(inputCoefficients.bottom.modeCoefficients, inputCoefficients.top.modeCoefficients)
end

function numHarmonics(inputCoefficients::InputCoefficients)
    return numHarmonics(inputCoefficients.bottom)
end


function inputFields2InputCoefficients( inputFields::InputFields, W₀::ElectricEigenvectors)::InputCoefficients
    bottomCoefficientSet = FieldSetStack2CoefficientSet(inputFields.bottom, W₀)
    topCoefficientSet = FieldSetStack2CoefficientSet(inputFields.top, W₀)
    return InputCoefficients(bottomCoefficientSet, topCoefficientSet)
end

function inputCoefficients2InputFields( inputCoefficients::InputCoefficients, W₀::ElectricEigenvectors)::OutputFields
    bottomFieldSet = modeCoefficientSet2FieldSet(inputCoefficients.bottom, W₀)
    topFieldSet = modeCoefficientSet2FieldSet(inputCoefficients.top, W₀)
    return OutputFields(bottomCoefficientSet, topCoefficientSet)
end

function outputFields2OutputCoefficients( outputFields::OutputFields, W₀::ElectricEigenvectors)::OutputCoefficients
    bottomCoefficientSet = FieldSetStack2CoefficientSet(outputFields.bottom, W₀)
    topCoefficientSet = FieldSetStack2CoefficientSet(outputFields.top, W₀)
    return OutputCoefficients(bottomCoefficientSet, topCoefficientSet)
end

function outputCoefficients2OutputFields( outputCoefficients::OutputCoefficients, W₀::ElectricEigenvectors)::OutputFields
    bottomFieldSet = modeCoefficientSet2FieldSet(outputCoefficients.bottom, W₀)
    topFieldSet = modeCoefficientSet2FieldSet(outputCoefficients.top, W₀)
    return OutputFields(bottomFieldSet, topFieldSet)
end



# TODO: Should change to E*conj(E) ?
# Calculate E² for the Array (nHarmonics X 3) of E fields
function calcE²(fields)
    @assert size(fields,2) == 3

    return abs.(fields[:,X]).^2 + abs.(fields[:,Y]).^2 + abs.(fields[:,Z]).^2
end
