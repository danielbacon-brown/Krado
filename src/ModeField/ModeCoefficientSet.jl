# Represents the mode coefficients travelling in a particular direction in a particular layer
mutable struct ModeCoefficientSet
    
    # Vector has length equal to the number of Harmonics
    # Values 1:nHarmonics are x-component coefficients
    # Values (nHarmonics+1):(2*nHarmonics) are y-component coefficients
    modeCoefficients::Vector{ComplexF64}
    
    # Defines whether this set represents modes travelling in the forward direction (+Z)
    isForward::Bool
    
    function ModeCoefficientSet(modeCoefficients::Vector{ComplexF64}, isForward::Bool)
        return new(modeCoefficients, isForward)
    end
    
end

function numHarmonics(modeCoefficientSet::ModeCoefficientSet)
    return half(length(modeCoefficientSet.modeCoefficients))
end
