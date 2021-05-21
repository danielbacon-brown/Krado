# Represents the relation between the mode coefficients and the mode fields 
# Denoted as W in the reference

# FROM REFERENCE:
# W – Square matrix who’s column
# vectors describe the “modes” that
# can exist in the material. These are
# essentially pictures of the modes
# which quantify the relative
# amplitudes of Ex, Ey, Hx, and Hy.



mutable struct ElectricEigenvectors

    # Array that is used to convert between mode fields and coefficients
    eigenvectors::Array{Float64,2}

    function ElectricEigenvectors(eigenvectors::Array{Float64,2})
        return new(eigenvectors)
    end
end

# Convert mode fields to mode coefficients
function FieldSetStack2CoefficientSet(FieldSetStack::FieldSetStack, W::ElectricEigenvectors)::ModeCoefficientSet
    coefficients = inv(W.eigenvectors) * FieldSetStack.modeFields
    return ModeCoefficientSet(coefficients, FieldSetStack.isForward)
end

# Convert mode coefficients to mode fields
function modeCoefficientSet2FieldSet(modeCoefficientSet::ModeCoefficientSet, W::ElectricEigenvectors)
    fields = W.eigenvectors * modeCoefficientSet.modeCoefficients
    return FieldSetStack(fields, modeCoefficientSet.isForward)
end
