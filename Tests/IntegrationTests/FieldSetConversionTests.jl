# module FieldSetConversions
using Test
include("../../src/IncludeKrado.jl")

@testset "FieldSetConversions" begin

println("starting")

# Define Materials and layers
matCol = MaterialCollection()
addMaterial!(matCol,"Air", Material(ConstantPermittivity(1)) )
layerBottom = SemiInfiniteLayerDefinition("Air")
layerTop = SemiInfiniteLayerDefinition("Air")
layerStack = [layerBottom, layerTop]

nAir = 1

lattice = Lattice([1.0,0.0], [0.0,1.0])

# Harmonics:
M,N = 0,0
harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)
# Results to calculate:
analysisDefinition = AllModesAnalysisDefinition()


# Define input boundary
wavenumber = WavenumberByλ₀(0.532*μm)
k₀ = getk₀(wavenumber)
θ = 1e-6
ϕ = 45*degrees
# A = [ 1, 1im ]
A = [ 1, 0 ]
B = [ 1, 1im ]
mainHarmonicOrder = [0,0]
isTop = false
AbyϖbottomA = Dict{_2VectorInt,_2VectorComplex}()
AbyϖbottomB = Dict{_2VectorInt,_2VectorComplex}()
Abyϖtop = Dict{_2VectorInt,_2VectorComplex}()
AbyϖbottomA[_2VectorInt(0,0)] = A
AbyϖbottomB[_2VectorInt(0,0)] = B
boundaryDefinitionA = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, mainHarmonicOrder, isTop, AbyϖbottomA, Abyϖtop)
boundaryDefinitionB = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, mainHarmonicOrder, isTop, AbyϖbottomB, Abyϖtop)

simulationDefinitionA = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinitionA, matCol, analysisDefinition)
simulationDefinitionB = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinitionB, matCol, analysisDefinition)

derivedParametersA = DerivedParameters(simulationDefinitionA)
derivedParametersB = DerivedParameters(simulationDefinitionB)
kVectorSet = derivedParametersA.kVectorSet

# K-vector:
kXYnorm = kVectorSet.kᵢNorm[1]
isForward = true
@test isapprox( kXYnorm, [1/sqrt(2), 0], rtol=1e-3, atol=1e-5)
kXYZ = kXYtokXYZ(kXYnorm*k₀, nAir, wavenumber, isForward)
@test isapprox( kXYZ/k₀, [1/sqrt(2), 0, 1/sqrt(2)], rtol=1e-3, atol=1e-5 )

# INPUT FIELDS:
# Create fields in terms of SP
fieldSetSPbottomA, fieldSetSPtopA = calcSPinputFields(derivedParametersA.boundaryConditions, derivedParametersA.harmonicsSet)
fieldSetSPbottomB, fieldSetSPtopB = calcSPinputFields(derivedParametersB.boundaryConditions, derivedParametersB.harmonicsSet)
# @test fieldsSPbottom.fields[1,:] ≈ [1, 1im]
@test fieldSetSPbottomA.fields[1,:] ≈ [1, 0]
@test fieldSetSPtopA.fields[1,:] ≈ [0, 0]
@test fieldSetSPbottomA.isForward == true
@test fieldSetSPtopA.isForward == false
@test fieldSetSPbottomB.fields[1,:] ≈ [1, 1im]
@test fieldSetSPtopB.fields[1,:] ≈ [0, 0]
@test fieldSetSPbottomB.isForward == true
@test fieldSetSPtopB.isForward == false

# Convert SP fields to xyz
fieldSetXYZbottomA, fieldSetXYZtopA = convertSPinputFieldsToXYZ( fieldSetSPbottomA, fieldSetSPtopA, kVectorSet, layerStack, matCol, wavenumber)
fieldSetXYZbottomB, fieldSetXYZtopB = convertSPinputFieldsToXYZ( fieldSetSPbottomB, fieldSetSPtopB, kVectorSet, layerStack, matCol, wavenumber)
@test isapprox( fieldSetXYZbottomA.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
@test fieldSetXYZbottomA.isForward == true
@test isapprox( fieldSetXYZbottomB.fields[1,:], [1/sqrt(2)*1im, -1, -1/sqrt(2)*1im], rtol=1e-3, atol=1e-5)
@test fieldSetXYZbottomB.isForward == true

# fig, ax = plot3Dlattice(simulationDefinitionB.lattice, simulationDefinitionB.layerStack; scale=1.0, title = "fieldSetXYZbottomB 1")
# fig = figure(title, figsize=(5,5))
# ax = Axes3D(fig)
# @show fieldSetXYZbottomB.fields[1,:]
@test isapprox( (kXYZ ⋅ fieldSetXYZbottomB.fields[1,:]), 0, rtol=1e-3, atol=1e-5)
# add3DKandPVectorToPlot(kXYZ, fieldSetXYZbottomB.fields[1,:], [0,0,0], wavenumber; scale=μm, Escale=1)

# Convert xyz fields to stacked fields
inputFieldStackBottomA = convertFieldSetXYZtoStack(fieldSetXYZbottomA)
inputFieldStackBottomB = convertFieldSetXYZtoStack(fieldSetXYZbottomB)
inputFieldStackTopA = convertFieldSetXYZtoStack(fieldSetXYZtopA)
inputFieldStackTopB = convertFieldSetXYZtoStack(fieldSetXYZtopB)
@test isapprox( inputFieldStackBottomA.modeFields, [0, -1], rtol=1e-3, atol=1e-5)
@test isapprox( inputFieldStackBottomB.modeFields, [1/sqrt(2)*1im, -1], rtol=1e-3, atol=1e-5)

# Convert stacked fields back to xyz fields
# kzNormBottom = derivedParametersA.kzNormBottom
fieldSetXYZbottomA = convertFieldSetStackToXYZ(inputFieldStackBottomA, kVectorSet, -derivedParametersA.kzNormBottom)
fieldSetXYZbottomB = convertFieldSetStackToXYZ(inputFieldStackBottomB, kVectorSet, -derivedParametersB.kzNormBottom)
# @show fieldsXYZbottom
@test isapprox( fieldSetXYZbottomA.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
@test fieldSetXYZbottomA.isForward == true
# @show kXYZ/k₀
# @show fieldSetXYZbottomB.fields[1,:]
@test isapprox( (kXYZ ⋅ fieldSetXYZbottomB.fields[1,:]), 0, rtol=1e-3, atol=1e-5)
# fig, ax = plot3Dlattice(simulationDefinitionB.lattice, simulationDefinitionB.layerStack; scale=1.0, title = "fieldSetXYZbottomB 2")
# add3DKandPVectorToPlot(kXYZ, fieldSetXYZbottomB.fields[1,:], [0,0,0], wavenumber; scale=μm, Escale=1)

@test isapprox( fieldSetXYZbottomB.fields[1,:], [1/sqrt(2)*1im, -1, -1/sqrt(2)*1im], rtol=1e-3, atol=1e-5)
@test fieldSetXYZbottomB.isForward == true

# convert xyz fields to SP fields
fieldsSPbottomA = convertFieldSetXYZtoSP( fieldSetXYZbottomA, kVectorSet, getMaterial(matCol,"Air"), wavenumber )
fieldsSPbottomB = convertFieldSetXYZtoSP( fieldSetXYZbottomB, kVectorSet, getMaterial(matCol,"Air"), wavenumber )
@test fieldsSPbottomA.fields[1,:] ≈ [1, 0]
@test fieldsSPbottomA.isForward == true
@test fieldsSPbottomB.fields[1,:] ≈ [1, 1im]
@test fieldsSPbottomB.isForward == true


# What is done in AllModesAnalysis:
inputFieldsA = InputFields(inputFieldStackBottomA, inputFieldStackTopA)
inputFieldsB = InputFields(inputFieldStackBottomB, inputFieldStackTopB)

# Done in vector plotting:
# Convert stacked fields back to xyz fields
inputFieldSetXYZbottomA = convertFieldSetStackToXYZ(inputFieldsA.bottom, kVectorSet, -derivedParametersA.kzNormBottom)
inputFieldSetXYZbottomB = convertFieldSetStackToXYZ(inputFieldsB.bottom, kVectorSet, -derivedParametersB.kzNormBottom)
@test isapprox( inputFieldSetXYZbottomA.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
@test inputFieldSetXYZbottomA.isForward == true
@test isapprox( inputFieldSetXYZbottomB.fields[1,:], [1/sqrt(2)*1im, -1, -1/sqrt(2)*1im], rtol=1e-3, atol=1e-5)
@test inputFieldSetXYZbottomB.isForward == true
@test isapprox( (kXYZ ⋅ inputFieldSetXYZbottomA.fields[1,:]), 0, rtol=1e-3, atol=1e-5)
@test isapprox( (kXYZ ⋅ inputFieldSetXYZbottomB.fields[1,:]), 0, rtol=1e-3, atol=1e-5)




### AllModesAnalysisDefinition:
# Calc common parameters
# derivedParametersA = DerivedParameters(simulationDefinitionA)
# # Calc input fields
# inputFields = calcInputFields(simulationDefinition, derivedParameters)
# @test isapprox( fieldsXYZbottom.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
# Calc global scattering matrix
SglobalA = calcGlobalScatteringMatrix(simulationDefinitionA, derivedParametersA)
SglobalB = calcGlobalScatteringMatrix(simulationDefinitionB, derivedParametersB)
# Propagate fields
# outputFields = propagateFields( Sglobal, inputFields, derivedParameters )
# Output data as a named tuple
# data = (inputFields = inputFields, outputFields = outputFields)

# Calculating output fields
outputFieldsA = propagateFields( SglobalA, inputFieldsA, derivedParametersA )
outputFieldsB = propagateFields( SglobalB, inputFieldsB, derivedParametersB )
@test isapprox( outputFieldsA.top.modeFields, inputFieldsA.bottom.modeFields, rtol=1e-3, atol=1e-5)
@test isapprox( outputFieldsB.top.modeFields, inputFieldsB.bottom.modeFields, rtol=1e-3, atol=1e-5)
# Convert output fields to XYZ
# Is negative by default, so must use a negative value for kz when calculating for bottom input.  Must also use a negative value for top output.
# @show derivedParametersA.kzNormBottom
outputFieldSetXYZbottomA = convertFieldSetStackToXYZ(outputFieldsA.bottom, kVectorSet, -derivedParametersA.kzNormBottom)
outputFieldSetXYZbottomB = convertFieldSetStackToXYZ(outputFieldsB.bottom, kVectorSet, -derivedParametersB.kzNormBottom)
outputFieldSetXYZtopA = convertFieldSetStackToXYZ(outputFieldsA.top, kVectorSet, derivedParametersA.kzNormTop )
outputFieldSetXYZtopB = convertFieldSetStackToXYZ(outputFieldsB.top, kVectorSet, derivedParametersB.kzNormTop )
@test isapprox( outputFieldSetXYZtopA.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
@test outputFieldSetXYZtopA.isForward == true
@test isapprox( outputFieldSetXYZtopB.fields[1,:], [1/sqrt(2)*1im, -1, -1/sqrt(2)*1im], rtol=1e-3, atol=1e-5)
@test outputFieldSetXYZtopA.isForward == true
@test isapprox( (kXYZ ⋅ outputFieldSetXYZtopA.fields[1,:]), 0, rtol=1e-3, atol=1e-5)
@test isapprox( (kXYZ ⋅ outputFieldSetXYZtopB.fields[1,:]), 0, rtol=1e-3, atol=1e-5)

# Convert output XYZ to SP
fieldsSPtopA = convertFieldSetXYZtoSP( outputFieldSetXYZtopA, kVectorSet, getMaterial(matCol,"Air"), wavenumber )
fieldsSPtopB = convertFieldSetXYZtoSP( outputFieldSetXYZtopB, kVectorSet, getMaterial(matCol,"Air"), wavenumber )
@test fieldsSPtopA.fields[1,:] ≈ [1, 0]
@test fieldsSPtopA.isForward == true
@test fieldsSPtopB.fields[1,:] ≈ [1, 1im]
@test fieldsSPtopB.isForward == true




# Following runSimulation:
derivedParametersA = DerivedParameters(simulationDefinitionA)
derivedParametersB = DerivedParameters(simulationDefinitionB)
allModeDataA = runSimulation(simulationDefinitionA)
allModeDataB = runSimulation(simulationDefinitionB)
@test isapprox( outputFieldsA.top.modeFields, inputFieldsA.bottom.modeFields, rtol=1e-3, atol=1e-5)
@test isapprox( outputFieldsB.top.modeFields, inputFieldsB.bottom.modeFields, rtol=1e-3, atol=1e-5)


# bottomOrders = [_2VectorInt(0,0),]
# topOrders = [_2VectorInt(0,0),]
# kVectorSet = derivedParametersA.kVectorSet
# harmonicsSet = derivedParams.harmonicsSet
wavenumber = getWavenumber(derivedParametersA.boundaryConditions)
# bottomOrderIndices = [ getOrderIndex(harmonicsSet, ϖ) for ϖ in bottomOrders]
# topOrderIndices = [ getOrderIndex(harmonicsSet, ϖ) for ϖ in topOrders]

inputFieldSetXYZbottomA = convertFieldSetStackToXYZ(allModeDataA.inputFields.bottom, kVectorSet, -derivedParametersA.kzNormBottom)
inputFieldSetXYZtopA = convertFieldSetStackToXYZ(allModeDataA.inputFields.top, kVectorSet, derivedParametersA.kzNormTop)
outputFieldSetXYZbottomA = convertFieldSetStackToXYZ(allModeDataA.outputFields.bottom, kVectorSet, derivedParametersA.kzNormBottom)
outputFieldSetXYZtopA = convertFieldSetStackToXYZ(allModeDataA.outputFields.top, kVectorSet, derivedParametersA.kzNormTop)
@test isapprox( outputFieldSetXYZtopA.fields, inputFieldSetXYZbottomA.fields, rtol=1e-3, atol=1e-5)
@test isapprox( outputFieldSetXYZtopA.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
@test isapprox( inputFieldSetXYZbottomA.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)

inputFieldSetXYZbottomB = convertFieldSetStackToXYZ(allModeDataB.inputFields.bottom, kVectorSet, -derivedParametersB.kzNormBottom)
inputFieldSetXYZtopB = convertFieldSetStackToXYZ(allModeDataB.inputFields.top, kVectorSet, derivedParametersB.kzNormTop)
outputFieldSetXYZbottomB = convertFieldSetStackToXYZ(allModeDataB.outputFields.bottom, kVectorSet, derivedParametersB.kzNormBottom)
outputFieldSetXYZtopB = convertFieldSetStackToXYZ(allModeDataB.outputFields.top, kVectorSet, derivedParametersB.kzNormTop)
@test isapprox( outputFieldSetXYZtopB.fields, inputFieldSetXYZbottomB.fields, rtol=1e-3, atol=1e-5)
@test isapprox( outputFieldSetXYZtopB.fields[1,:], [1/sqrt(2)*1im, -1, -1/sqrt(2)*1im], rtol=1e-3, atol=1e-5)
@test isapprox( inputFieldSetXYZbottomB.fields[1,:], [1/sqrt(2)*1im, -1, -1/sqrt(2)*1im], rtol=1e-3, atol=1e-5)

inputFieldSetSPbottomA = convertFieldSetXYZtoSP( inputFieldSetXYZbottomA, kVectorSet, getMaterial(matCol,"Air"), wavenumber )
outputFieldSetSPtopA = convertFieldSetXYZtoSP( outputFieldSetXYZtopA, kVectorSet, getMaterial(matCol,"Air"), wavenumber )
@test inputFieldSetSPbottomA.fields[1,:] ≈ [1, 0]
@test inputFieldSetSPbottomA.isForward == true
@test outputFieldSetSPtopA.fields[1,:] ≈ [1, 0]
@test outputFieldSetSPtopA.isForward == true

inputFieldSetSPbottomB = convertFieldSetXYZtoSP( inputFieldSetXYZbottomB, kVectorSet, getMaterial(matCol,"Air"), wavenumber )
outputFieldSetSPtopB = convertFieldSetXYZtoSP( outputFieldSetXYZtopB, kVectorSet, getMaterial(matCol,"Air"), wavenumber )
@test inputFieldSetSPbottomB.fields[1,:] ≈ [1, 1im]
@test inputFieldSetSPbottomB.isForward == true
@test outputFieldSetSPtopB.fields[1,:] ≈ [1, 1im]
@test outputFieldSetSPtopB.isForward == true

inputFieldSetXYZbottomA, inputFieldSetXYZtopA, inputFieldSetSPbottomA, inputFieldSetSPtopA = analyzeInputFields(inputFieldsA, derivedParametersA, simulationDefinitionA.layerStack, simulationDefinitionA.materialCollection, getWavenumber(simulationDefinitionA) )
inputFieldSetXYZbottomB, inputFieldSetXYZtopB, inputFieldSetSPbottomB, inputFieldSetSPtopB = analyzeInputFields(inputFieldsB, derivedParametersB, simulationDefinitionB.layerStack, simulationDefinitionB.materialCollection, getWavenumber(simulationDefinitionB) )
@test isapprox( inputFieldSetXYZbottomA.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
@test inputFieldSetXYZbottomA.isForward == true
@test isapprox( inputFieldSetXYZbottomB.fields[1,:], [1/sqrt(2)*1im, -1, -1/sqrt(2)*1im], rtol=1e-3, atol=1e-5)
@test inputFieldSetXYZbottomB.isForward == true
@test isapprox( (kXYZ ⋅ inputFieldSetXYZbottomA.fields[1,:]), 0, rtol=1e-3, atol=1e-5)
@test isapprox( (kXYZ ⋅ inputFieldSetXYZbottomB.fields[1,:]), 0, rtol=1e-3, atol=1e-5)

@test inputFieldSetSPbottomA.fields[1,:] ≈ [1, 0]
@test inputFieldSetSPbottomA.isForward == true
@test outputFieldSetSPtopA.fields[1,:] ≈ [1, 0]
@test outputFieldSetSPtopA.isForward == true
@test inputFieldSetSPbottomB.fields[1,:] ≈ [1, 1im]
@test inputFieldSetSPbottomB.isForward == true
@test outputFieldSetSPtopB.fields[1,:] ≈ [1, 1im]
@test outputFieldSetSPtopB.isForward == true




end;

# end # module
