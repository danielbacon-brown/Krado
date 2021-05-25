# Contains the 3-dimensional, XYZ fields values corresponding to modes travelling in a particular layer in a particular direction.
mutable struct FieldSetXYZ

    # Fields array is of size (# harmonics, 3)
    fields::Array{ComplexF64,2}

    isForward::Bool

    function FieldSetXYZ(fields::Array{ComplexF64,2}, isForward::Bool)
        return new(fields, isForward)
    end

end

# Calculates the XYZ field data based on the XY field data.
# function convertFieldSetStackToXYZ(fieldSetStack::FieldSetStack, kVectorSet::KVectorSet, kz::Vector{T}) where T<:Number
#
#     k₀ = getk₀(kVectorSet.wavenumber)
#     kzDirection = bool2posNeg(fieldSetStack.isForward) # Reverse value if backwards propagating
#
#     nHarmonics = numHarmonics(kVectorSet)
#     fieldsXYZ = Array{ComplexF64,2}(undef, (nHarmonics,3))
#
#     fieldsXYZ[:,X] = fieldSetStack.modeFields[1:nHarmonics]
#     fieldsXYZ[:,Y] = fieldSetStack.modeFields[ (nHarmonics+1):(2*nHarmonics)]
#     fieldsXYZ[:,Z] = -inv(Diagonal(kz))*(kVectorSet.KxNorm*fieldsXYZ[:,X] + kVectorSet.KyNorm*fieldsXYZ[:,Y])
#     return FieldSetXYZ(fieldsXYZ, fieldSetStack.isForward)
# end

function calcE²(fieldSetXYZ::FieldSetXYZ)
    return abs.(fieldSetXYZ.fields[:,X]).^2 + abs.(fieldSetXYZ.fields[:,Y]).^2 + abs.(fieldSetXYZ.fields[:,Z]).^2
end


function calcPowerFlux(fieldSetXYZ::FieldSetXYZ, kz::Vector{ComplexF64})
    powerFluxes = c*ϵ₀/2 * (kz .* calcE²(fieldSetXYZ))
    return powerFluxes
end





function getField3ByOrder(fieldSetXYZ::FieldSetXYZ, order::Int64)
    return _3VectorComplex(fieldSetXYZ.fields[order,:])
end

function numHarmonics(fieldSetXYZ::FieldSetXYZ)
    return size(fieldSetXYZ.fields, 1)
end
