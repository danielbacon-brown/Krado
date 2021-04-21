# Describes the light conditions of light entering and leaving the stack
abstract type BoundaryConditions end

# Input is defined by a single k-vector defining the 0-order mode, and a dictionary of pairs of vectors and amplitudes
mutable struct InputByOrderBoundaryConditions <: BoundaryConditions

    # Vacuum wavenumber
    wavenumber::Wavenumber

    # k 2-vector of the first 0-order mode
    kXY₀::_2VectorFloat

    # Dictionary defining the mode amplitudes, A, for each integer harmonic order, ϖ
    AbyϖBottom::Dict{_2VectorInt, _2VectorComplex} # side 1 = first layer
    AbyϖTop::Dict{_2VectorInt, _2VectorComplex} # side 2 = last layer

    function InputByOrderBoundaryConditions(wavenumber::Wavenumber, kXY₀::_2VectorFloat, AbyϖBottom::Dict{_2VectorInt,_2VectorComplex}, AbyϖTop::Dict{_2VectorInt,_2VectorComplex})
        new(wavenumber, kXY₀, AbyϖBottom, AbyϖTop)
    end

end

function InputByOrderBoundaryConditions( boundaryDef::InputByOrderBoundaryDefinition, matCol::MaterialCollection, layerStack)
    n = getn( getBoundaryLayer(layerStack, boundaryDef.isTop), matCol, boundaryDef.wavenumber)
    if boundaryDef.isTop
        kNormal = -1*normalIncidenceKvector(boundaryDef.wavenumber)
    else
        kNormal = normalIncidenceKvector(boundaryDef.wavenumber)
    end

    # Calculate z-order k-vector.
    kXYZ = rotateVector(kNormal, boundaryDef.θ, boundaryDef.ϕ) * n
    kXY₀ = _2VectorFloat(0,0)
    try
        kXY₀ = _2VectorFloat(getXY(kXYZ))
    catch InexactError
        @warn "k-vector has a complex component due to absorption in the material in which the light is injected.  If this absorption is high, the results may be innacurate.  Forcing this value to be real."
        kXY₀ = _2VectorFloat(real.(getXY(kXYZ)))
    end

    return InputByOrderBoundaryConditions(boundaryDef.wavenumber, kXY₀, boundaryDef.Abyϖbottom, boundaryDef.Abyϖtop)
end


function createKVectorSet(boundaryDefinition::InputByOrderBoundaryDefinition, boundaryConditions::InputByOrderBoundaryConditions, GvectorSet::GvectorSet)
    return createKVectorSet(boundaryDefinition.wavenumber, boundaryConditions.kXY₀, boundaryDefinition.mainHarmonicOrder, GvectorSet)
end

# Returns the 3-vector of the zero-order k-vector
function getkXYZ₀( boundaryConditions::InputByOrderBoundaryConditions, layerStack, matCol::MaterialCollection, kzPositive)
    wavenumber = getWavenumber(boundaryConditions)
    if ~kzPositive
        n = getn( last(layerStack), matCol, wavenumber)
    else
        n = getn( first(layerStack), matCol, wavenumber)
    end

    return kXYtokXYZ(boundaryConditions.kXY₀, n, wavenumber, kzPositive )
end


# Define a k-vector set using the k-vector of the boundary conditions
function createKVectorSet(Gvectors::GvectorSet, boundCond::InputByOrderBoundaryConditions)
    ϖ = _2VectorInt(0,0)
    return createKVectorSet(boundCond.wavenumber, boundCond.kXY₀, ϖ, Gvectors)
end


# Returns FieldSetSPs for the top layer (backward incidence) and bottom layer (forward incidence) based corresponding to the boundary conditions.
function calcSPinputFields(boundaryConditions::InputByOrderBoundaryConditions, harmonicsSet::HarmonicsSet)

    modeFieldSPbottom = calcFieldSetSP( boundaryConditions.AbyϖBottom, harmonicsSet, FORWARD)
    modeFieldSPtop = calcFieldSetSP( boundaryConditions.AbyϖTop, harmonicsSet, BACKWARD)

    return modeFieldSPbottom, modeFieldSPtop
end


# Converts the bottom and top SP field sets to XYZ field sets.  Requires calculation of n of top and bottom layers.
function convertSPinputFieldsToXYZ(modeFieldsSPbottom::FieldSetSP, modeFieldsSPtop::FieldSetSP, kVectorSet::KVectorSet, layerStack, matCol::MaterialCollection, wavenumber::Wavenumber)

    bottomMaterial = getMaterial(matCol, first(layerStack).backgroundMaterialName)
    modeFieldsXYZbottom = convertFieldSetSPtoXYZ( modeFieldsSPbottom, kVectorSet, bottomMaterial, wavenumber )
    topMaterial = getMaterial(matCol, last(layerStack).backgroundMaterialName)
    modeFieldsXYZtop = convertFieldSetSPtoXYZ( modeFieldsSPtop, kVectorSet, topMaterial, wavenumber )
    return modeFieldsXYZbottom, modeFieldsXYZtop
end

# Converts the bottom and top SP field sets to XYZ field sets.  Requires calculation of n of top and bottom layers.
function convertXYZoutputFieldsToSP(fieldSetXYZbottom::FieldSetXYZ, fieldSetXYZtop::FieldSetXYZ, kVectorSet::KVectorSet, layerStack, matCol::MaterialCollection, wavenumber::Wavenumber)

    bottomMaterial = getMaterial(matCol, first(layerStack).backgroundMaterialName)
    fieldSetSPbottom = convertFieldSetXYZtoSP( fieldSetXYZbottom, kVectorSet, bottomMaterial, wavenumber )
    topMaterial = getMaterial(matCol, last(layerStack).backgroundMaterialName)
    fieldSetSPtop = convertFieldSetXYZtoSP( fieldSetXYZtop, kVectorSet, topMaterial, wavenumber )
    return fieldSetSPbottom, fieldSetSPtop
end




# Returns InputFields representing stacked field sets for incidence from bottom and top.
function calcInputFields(boundaryConditions::InputByOrderBoundaryConditions, harmonicsSet::HarmonicsSet, kVectorSet::KVectorSet, layerStack, matCol::MaterialCollection, wavenumber::Wavenumber )::InputFields

    # Create fields in terms of SP
    modeFieldsSPbottom, modeFieldsSPtop = calcSPinputFields(boundaryConditions, harmonicsSet)

    # Convert SP fields to xyz
    modeFieldsXYZbottom, modeFieldsXYZtop = convertSPinputFieldsToXYZ(modeFieldsSPbottom, modeFieldsSPtop, kVectorSet, layerStack, matCol, wavenumber)

    # Convert xyz fields to stacked fields
    inputFieldsBottom = convertFieldSetXYZtoStack(modeFieldsXYZbottom)
    inputFieldsTop = convertFieldSetXYZtoStack(modeFieldsXYZtop)

    return InputFields(inputFieldsBottom, inputFieldsTop)
end





function BoundaryConditions( boundaryDef::InputByOrderBoundaryDefinition, simulationDef::SimulationDefinition)
    return InputByOrderBoundaryConditions( boundaryDef, simulationDef.materialCollection, simulationDef.layerStack)
end


function getWavenumber(inputByOrderBoundaryConditions::InputByOrderBoundaryConditions)
    return inputByOrderBoundaryConditions.wavenumber
end
