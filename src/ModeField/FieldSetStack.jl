# Represents the x,y fields for modes in a particular direction in a particular layer
mutable struct FieldSetStack

    # Vector has length equal to the number of Harmonics
    # Values 1:nHarmonics are x-component coefficients
    # Values (nHarmonics+1):(2*nHarmonics) are y-component coefficients
    modeFields::Vector{ComplexF64}

    # Defines whether this set represents modes travelling in the forward direction (+Z)
    isForward::Bool
    
    function FieldSetStack(modeCoefficients::Vector{ComplexF64}, isForward::Bool)
        return new(modeCoefficients, isForward)
    end

end

function numHarmonics(fieldSetStack::FieldSetStack)
    return half( size(fieldSetStack.fields, 1) )
end
