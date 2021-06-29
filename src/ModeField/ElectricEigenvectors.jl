# Represents the relation between the mode coefficients and the mode fields 
# Denoted as W in the reference
mutable struct ElectricEigenvectors{T} 
    
    # Array that is used to convert between mode fields and coefficients
    matrix::Array{Complex{T},2}
    
    function ElectricEigenvectors{T}(matrix::Array{W,2}) where {T<:Real, W<:Number}
        return new( convert(Array{Complex{T},2}, matrix) )
    end
end

function Base.convert( ::Type{ElectricEigenvectors{T1}}, eigenvectors::ElectricEigenvectors{T2}) where {T1<:Real, T2<:Real}
    return ElectricEigenvectors{T1}( convert(Array{Complex{T1},2}, eigenvectors.matrix) )
end


# Convert mode fields to mode coefficients
function FieldSetStack2CoefficientSet(fieldSetStack::FieldSetStack, W::ElectricEigenvectors)::ModeCoefficientSet
    coefficients = inv(W.matrix) * fieldSetStack.modeFields
    return ModeCoefficientSet(coefficients, fieldSetStack.isForward) 
end

# Convert mode coefficients to mode fields
function modeCoefficientSet2FieldSet(modeCoefficientSet::ModeCoefficientSet, W::ElectricEigenvectors)
    fields = W.matrix * modeCoefficientSet.modeCoefficients    
    return FieldSetStack(fields, modeCoefficientSet.isForward)
end


# Denoted as V in the reference
mutable struct MagneticEigenvectors{T}
    
    # Array that is used to convert between mode fields and coefficients
    matrix::Array{Complex{T},2}
    
    function MagneticEigenvectors{T}(matrix::Array{W,2}) where {T<:Real, W<:Number}
        return new( convert(Array{Complex{T},2}, matrix) )
    end
end

function Base.convert( ::Type{MagneticEigenvectors{T1}}, eigenvectors::MagneticEigenvectors{T2}) where {T1<:Real, T2<:Real}
    return MagneticEigenvectors{T1}( convert(Array{Complex{T1},2}, eigenvectors.matrix) )
end
