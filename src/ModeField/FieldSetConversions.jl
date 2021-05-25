# FieldSetConversions

function convertFieldSetStackToXYZ(fieldSetStack::FieldSetStack, kVectorSet::KVectorSet, kz::Vector{T}) where T<:Number

    k₀ = getk₀(kVectorSet.wavenumber)
    kzDirection = bool2posNeg(fieldSetStack.isForward) # Reverse value if backwards propagating

    nHarmonics = numHarmonics(kVectorSet)
    fieldsXYZ = Array{ComplexF64,2}(undef, (nHarmonics,3))

    fieldsXYZ[:,X] = fieldSetStack.modeFields[1:nHarmonics]
    fieldsXYZ[:,Y] = fieldSetStack.modeFields[ (nHarmonics+1):(2*nHarmonics)]
    #UNSURE
    fieldsXYZ[:,Z] = -inv(Diagonal(kz))*(kVectorSet.KxNorm*fieldsXYZ[:,X] + kVectorSet.KyNorm*fieldsXYZ[:,Y]) # old.  Problem with imaginary fields and p-polarizaiton
    # fieldsXYZ[:,Z] = -inv.(Diagonal(kz))*(kVectorSet.KxNorm*fieldsXYZ[:,X] + kVectorSet.KyNorm*fieldsXYZ[:,Y])
    # fieldsXYZ[:,Z] = -conj( inv(Diagonal(kz))*(kVectorSet.KxNorm*fieldsXYZ[:,X] + kVectorSet.KyNorm*fieldsXYZ[:,Y]) ) # new.  Fixes the conversion from Stack to XYZ
    # fieldsXYZ[:,Z] = -inv(conj(Diagonal(kz)))*(kVectorSet.KxNorm*fieldsXYZ[:,X] + kVectorSet.KyNorm*fieldsXYZ[:,Y])
    return FieldSetXYZ(fieldsXYZ, fieldSetStack.isForward)
end

# Converts the fieldSet in XYZ to XY with a stacked format.  Allows it to be converted to mode coefficients.
function convertFieldSetXYZtoStack(fieldSetXYZ::FieldSetXYZ )

    nHarmonics = numHarmonics(fieldSetXYZ)

    fieldStack = Vector{ComplexF64}(undef, nHarmonics*2)

    for ϖindex = 1:nHarmonics
        Pxy = getXY(fieldSetXYZ.fields[ϖindex,:])
        fieldStack[ϖindex] = Pxy[X]
        fieldStack[nHarmonics+ϖindex] = Pxy[Y]
    end

    return FieldSetStack(fieldStack, fieldSetXYZ.isForward)
end


# Converts the given SP field set to an XYZ field set assuming that the fields are in the given material with the given frequency.
function convertFieldSetSPtoXYZ( fieldSetSP::FieldSetSP, kVectorSet::KVectorSet, material::AbstractMaterial, wavenumber::Wavenumber )

    nHarmonics = numHarmonics(kVectorSet)
    isForward = fieldSetSP.isForward

    n = convert_ϵ2n(calc_ϵ( material, wavenumber))

    fieldsXYZ = Array{ComplexF64,2}(undef,(nHarmonics,3))

    for ϖindex in 1:nHarmonics
        kXY = getkXYnorm(kVectorSet, ϖindex) * getk₀(kVectorSet.wavenumber)
        kXYZ = kXYtokXYZ(kXY, n, wavenumber, isForward)
        fieldSP = fieldSetSP.fields[ϖindex,:]
        fieldsXYZ[ϖindex,:] = fieldSPtoFieldXYZ(kXYZ, fieldSP)
    end

    return FieldSetXYZ(fieldsXYZ, isForward)
end


# Converts the given XYZ field set to an SP field set assuming that the fields are in the given material with the given frequency.
function convertFieldSetXYZtoSP( fieldSetXYZ::FieldSetXYZ, kVectorSet::KVectorSet, material::AbstractMaterial, wavenumber::Wavenumber )

    nHarmonics = numHarmonics(kVectorSet)
    isForward = fieldSetXYZ.isForward

    n = convert_ϵ2n(calc_ϵ( material, wavenumber))

    fieldsSP = Array{ComplexF64,2}(undef,(nHarmonics,2))

    for ϖindex in 1:nHarmonics
        kXYnorm = getkXYnorm(kVectorSet, ϖindex)
        kXY = kXYnorm * getk₀(wavenumber)
        kXYZ = kXYtokXYZ(kXY, n, wavenumber, isForward)
        fieldXYZ = fieldSetXYZ.fields[ϖindex,:]
        fieldsSP[ϖindex,:] = fieldXYZtoFieldSP(kXYZ, fieldXYZ)
    end

    return FieldSetSP(fieldsSP, isForward)
end
