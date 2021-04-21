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



# TODO: GET RID OF THIS?
# Contains data for both the forward and backward. For X,Y,Z.
# Forward fields correspond to modes travelling in the +Z direction.
# Backward fields correspond to modes travelling in the -Z direction
mutable struct LayerFields
    
    forward::Array{ComplexF64,2}
    backward::Array{ComplexF64,2}
    
    function LayerFields(forward::Array{ComplexF64,2}, backward::Array{ComplexF64,2})
        new(forward, backward)
    end
end

# TODO: REPLACE INTIRELY WITH FIELDSETS 
# function stackedFieldsData2joinedFieldsData(fieldsXY::Array{ComplexF64,2}, kz::T1, kVectorSet::KVectorSet) where T1<:Union{ Array{ComplexF64,2}, LinearAlgebra.Diagonal{Complex{Float64},Array{Complex{Float64},1}} }
# 
#     @assert size(fieldsXY,2) == 1
#     nHarmonics = numHarmonics(kVectorSet)
#     fields = Array{ComplexF64,2}(undef, (nHarmonics,3))
# 
#     fields[:,X] = fieldsXY[1:nHarmonics]
#     fields[:,Y] = fieldsXY[(nHarmonics+1):(2*nHarmonics)]
#     fields[:,Z] = -1*inv(kz)*(kVectorSet.Kx*fields[:,X] + kVectorSet.Ky*fields[:,Y])
#     return fields
# end

# TODO: NOT USED
# Converts the input and output XY field data into XYZ data for each layer0 
# function getFieldsBottom(inputFields::InputFields, outputFields::OutputFields, kzᵦ::T1,  kVectorSet::KVectorSet)::LayerFields where T1<:Union{ Array{ComplexF64,2}, LinearAlgebra.Diagonal{ Complex{Float64},Array{Complex{Float64},1}} }
# 
#     forwardFields = stackedFieldsData2joinedFieldsData(inputFields.bottom, kzᵦ, kVectorSet)
#     backwardFields = stackedFieldsData2joinedFieldsData(outputFields.bottom, -1*kzᵦ, kVectorSet)
# 
#     return LayerFields(forwardFields, backwardFields)
# end

# NOT USED
# Converts the input and output XY field data into XYZ data for each layer0 
# function getFieldsTop(inputFields::InputFields, outputFields::OutputFields, kzₜ::T1,  kVectorSet::KVectorSet)::LayerFields where T1<:Union{ Array{ComplexF64,2}, LinearAlgebra.Diagonal{Complex{Float64},Array{Complex{Float64},1}} }
# 
#     nHarmonics = numHarmonics(kVectorSet)
# 
#     forwardFields = stackedFieldsData2joinedFieldsData(outputFields.top, kzₜ, kVectorSet)
#     backwardFields = stackedFieldsData2joinedFieldsData(inputFields.top, -1*kzₜ, kVectorSet)
# 
#     return LayerFields(forwardFields, backwardFields)
# end

# TODO: Should change to E*conj(E) ???
# Calculate E² for the Array (nHarmonics X 3) of E fields
function calcE²(fields)
    @assert size(fields,2) == 3
    
    return abs.(fields[:,X]).^2 + abs.(fields[:,Y]).^2 + abs.(fields[:,Z]).^2
end


# function calcPowerFluxes(inputFields::InputFields, outputFields::OutputFields, derivedParameters::DerivedParameters)
# 
#     # Convert to 3-vector data
#     bottomFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzBottom)
#     topFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.top, derivedParameters.kVectorSet, derivedParameters.kzTop)
#     bottomFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.bottom, derivedParameters.kVectorSet, derivedParameters.kzBottom)
#     topFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.top, derivedParameters.kVectorSet, derivedParameters.kzTop)
# 
#     # Get power flux of each
#     inputBottomPowerFlux = calcPowerFlux(bottomFieldsInput, derivedParameters.kzBottom)
#     outputBottomPowerFlux = calcPowerFlux(bottomFieldsOutput, derivedParameters.kzBottom)
#     inputTopPowerFlux = calcPowerFlux(topFieldsInput, derivedParameters.kzBottom)
#     outputTopPowerFlux = calcPowerFlux(topFieldsOutput, derivedParameters.kzTop)
# 
#     totalInputPowerFlux = sum(real(inputBottomPowerFlux)) + sum(real(inputTopPowerFlux))
# 
#     return inputBottomRelativeFlux, inputTopRelativeFlux, outputBottomRelativeFlux, outputTopRelativeFlux 
# end
