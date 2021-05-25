module FieldSetConversions
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
@show fieldSetXYZbottomB.fields[1,:]
@test (kXYZ ⋅ fieldSetXYZbottomB.fields[1,:]) ≈ 0
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
fieldSetXYZbottomA = convertFieldSetStackToXYZ(inputFieldStackBottomA, kVectorSet, derivedParametersA.kzNormBottom)
fieldSetXYZbottomB = convertFieldSetStackToXYZ(inputFieldStackBottomB, kVectorSet, derivedParametersB.kzNormBottom)
# @show fieldsXYZbottom
@test isapprox( fieldSetXYZbottomA.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
@test fieldSetXYZbottomA.isForward == true
# @show kXYZ/k₀
@show fieldSetXYZbottomB.fields[1,:]
@test (kXYZ ⋅ fieldSetXYZbottomB.fields[1,:]) ≈ 0
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
inputFieldSetXYZbottomA = convertFieldSetStackToXYZ(inputFieldsA.bottom, kVectorSet, derivedParametersA.kzNormBottom)
inputFieldSetXYZbottomB = convertFieldSetStackToXYZ(inputFieldsB.bottom, kVectorSet, derivedParametersB.kzNormBottom)
@test isapprox( inputFieldSetXYZbottomA.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
@test inputFieldSetXYZbottomA.isForward == true
@test isapprox( inputFieldSetXYZbottomB.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
@test inputFieldSetXYZbottomB.isForward == true





# Testing only relative values:
# wavenumber = WavenumberByλ₀(0.532*μm)
# k₀ = getk₀(wavenumber)
# θ = 1e-6
# ϕ = 45*degrees
# A = [ 1, 1im ]
# # A = [ 1, 0 ]
# mainHarmonicOrder = [0,0]
# isTop = false
# Abyϖbottom = Dict{_2VectorInt,_2VectorComplex}()
# Abyϖtop = Dict{_2VectorInt,_2VectorComplex}()
# Abyϖbottom[_2VectorInt(0,0)] = A
# boundaryDefinition = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, mainHarmonicOrder, isTop, Abyϖbottom, Abyϖtop)
#
# simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition)
#
# derivedParameters = DerivedParameters(simulationDefinition)
# kVectorSet = derivedParameters.kVectorSet
#
#
#
#
#
# # Done with add3DinjectedKandPVectorsToPlot
# fieldSetSPbottom, fieldSetSPtop = calcSPinputFields( derivedParameters.boundaryConditions, derivedParameters.harmonicsSet)
# fieldSetXYZbottom, fieldSetXYZtop = convertSPinputFieldsToXYZ( fieldSetSPbottom, fieldSetSPtop, derivedParameters.kVectorSet, simulationDefinition.layerStack, simulationDefinition.materialCollection, wavenumber)
# injectedOrderIndicesBottom = getInjectedOrderIndices(fieldSetXYZbottom)
# # add3DKandPVectorsToPlot( fieldSetXYZbottom, injectedOrderIndicesBottom, BOTTOM, simulationDef, derivedParams; scale=scale, Escale=Escale)
# injectedOrderIndicesTop = getInjectedOrderIndices(fieldSetXYZtop)
# # add3DKandPVectorsToPlot( fieldSetXYZtop, injectedOrderIndicesTop, TOP, simulationDef, derivedParams; scale=scale, Escale=Escale)
#
# # Done with add3DlistedKandPVectorsToPlot
# kVectorSet = derivedParameters.kVectorSet
# harmonicsSet = derivedParameters.harmonicsSet
# wavenumber = getWavenumber(derivedParameters.boundaryConditions)
# inputFieldSetXYZbottom = convertFieldSetStackToXYZ(inputFields.bottom, kVectorSet, derivedParameters.kzNormBottom)
# inputFieldSetXYZtop = convertFieldSetStackToXYZ(inputFields.top, kVectorSet, derivedParameters.kzNormTop)
# # outputFieldSetXYZbottom = convertFieldSetStackToXYZ(outputFields.bottom, kVectorSet, derivedParameters.kzNormBottom)
# # outputFieldSetXYZtop = convertFieldSetStackToXYZ(outputFields.top, kVectorSet, derivedParameters.kzNormTop)
# # bottomOrderIndices = [ getOrderIndex(harmonicsSet, ϖ) for ϖ in bottomOrders]
# # topOrderIndices = [ getOrderIndex(harmonicsSet, ϖ) for ϖ in topOrders]
# # add3DKandPVectorsToPlot( inputFieldSetXYZbottom, bottomOrderIndices, BOTTOM, simulationDef, derivedParams; scale=scale, Escale=Escale)
# # add3DKandPVectorsToPlot( inputFieldSetXYZtop, topOrderIndices, TOP, simulationDef, derivedParams; scale=scale, Escale=Escale)
#
# # Compare:
# @test isapprox( fieldSetXYZbottom.fields, inputFieldSetXYZbottom.fields,  rtol=1e-3, atol=1e-5)
# @test fieldSetXYZbottom.isForward == inputFieldSetXYZbottom.isForward





### AllModesAnalysisDefinition:
# Calc common parameters
derivedParameters = DerivedParameters(simulationDefinition)
# Calc input fields
inputFields = calcInputFields(simulationDefinition, derivedParameters)
@test isapprox( fieldsXYZbottom.fields[1,:], [0, -1, 0], rtol=1e-3, atol=1e-5)
# Calc global scattering matrix
Sglobal = calcGlobalScatteringMatrix(simulationDefinition, derivedParameters)
# Propagate fields
outputFields = propagateFields( Sglobal, inputFields, derivedParameters )
# Output data as a named tuple
# data = (inputFields = inputFields, outputFields = outputFields)




end;

end # module
