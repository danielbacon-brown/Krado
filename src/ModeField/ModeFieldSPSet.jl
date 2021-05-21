# Represents the SP fields for modes in a particular direction in a particular layer
mutable struct FieldSetSP

    # Array containing field data.  Size is (# harmonics, 2).  1st column is S.  2nd column is P.
    fields::Array{ComplexF64,2}

    # Defines whether this set represents modes travelling in the forward direction (+Z)
    isForward::Bool

    function FieldSetSP(fields::Array{<:Complex,2}, isForward::Bool)
        return new(fields, isForward)
    end

end



# Returns a FieldSetSP based corresponding to the Abyϖ.  Abyϖ is in terms of SP.
function calcFieldSetSP( Abyϖ::Dict{_2VectorInt, _2VectorComplex}, harmonicsSet::HarmonicsSet, isForward::Bool)::FieldSetSP

    fieldSetSP = zeros(ComplexF64, (numHarmonics(harmonicsSet),2))

    for (ϖ,A) in Abyϖ
        ϖindex = getOrderIndex( harmonicsSet, ϖ)
        fieldSetSP[ϖindex,:] = A
    end

    return FieldSetSP( fieldSetSP, isForward )
end

function numHarmonics(fieldSetSP::FieldSetSP)
    return size(fieldSetSP.fields, 1)
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
