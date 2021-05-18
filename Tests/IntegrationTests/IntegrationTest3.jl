module IntegrationTest3
# Using RCWA-Benchmark-Data-3x3 by Raymond Rumpf from empossible.net

using Test
# using LinearAlgebra
# using StaticArrays
# using Printf


println()
println()
println("Start:")

@testset "IntegrationTest3" begin
include("IntegrationTest3BenchmarkData.jl")

println()
println("Beginning IntegrationTest3")

include("../../src/IncludeKrado.jl")





######## DEFINE BOUNDARYDEFINITION #################################################################
# Incident wavevector
λ₀ = 2*cm
wavenumber = WavenumberByλ₀(λ₀)

# Benchmark appears to use nonstandard rotation method.  Here, θ is polar angle (rotation around z-axis) and ϕ is zenith angle (rotation around Y-axis).  ϕ rotation occurs first.
θ = 30 * degrees
ϕ = 60 * degrees
Eₚ = 0.70711
Eₛ = -0.70711im
inputAmplitudes = _2VectorComplex(Eₚ, Eₛ)
boundaryDefinition = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)


######## DEFINE MATERIALCOLLECTION #################################################################
matCol = MaterialCollection()

# Reflection material:
bottomPermittivity = 2.0
bottomPermeability = 1.0
addMaterial!(matCol,"bottom material", Material( ConstantPermittivity(bottomPermittivity), ConstantPermeability(bottomPermeability) ) )

# Transmission material:
topPermittivity = 9.0
topPermeability = 1.0
addMaterial!(matCol,"top material", Material( ConstantPermittivity(topPermittivity), ConstantPermeability(topPermeability) ) )

# Device material:
devicePermittivity = 6.0
devicePermeability = 1.0
addMaterial!(matCol,"device material", Material( ConstantPermittivity(devicePermittivity), ConstantPermeability(devicePermeability) ) )

######## DEFINE LATTICE ############################################################################
Lx = 1.75 * cm
Ly = 1.5 * cm
U̅ = [Lx, 0]
V̅ = [0, Ly]
lattice = Lattice(U̅, V̅)

######## DEFINE LAYERS #############################################################################
# Define layers.  Layer 1 contains triangle
layer1depth = 0.5 * cm  # thickness of layer 1
layer2depth = 0.3 * cm

# Define unpatterned layer:
layer2 = UniformLayerDefinition(layer2depth, "device material")

# Define shape: equilateral triangle pointed toward +y.
# Define a grid of materialNames so that we can define a GridLayerPattern
w = 0.8*Ly  # side-length of triangle
Nx = 512
Ny = round(Int64,Nx*Ly/Lx) # Number of pixels in y-direction

dx = Lx/Nx
dy = Ly/Ny # Y-length of each pixel

materialNameGrid = String["device material" for ix=1:Nx, iy=1:Ny ]

h = 0.5*sqrt(3)*w # height of triangle
ny = round(Int64,h/dy) # Number of rows for triangle
ny1 = round(Int64,(Ny - ny)/2) # top y point of triangle
ny2 = ny1 + ny - 1 # bottom y edge of triangle
for ny = ny1:ny2
    f = (ny-ny1)/(ny2-ny1)
    nx = round(Int64,f*w/Lx*Nx)
    nx1 = 1 + floor(Int64,(Nx-nx)/2)
    nx2 = nx1 + nx
    materialNameGrid[nx1:nx2, ny] = String["bottom material" for ix=nx1:nx2]
end
layerPattern1 = GridLayerPattern( materialNameGrid )

# FFT real-space coordinates
Nx = 512
Ny = round(Int64, Nx*Ly/Lx)
layerDivisions = [Nx, Ny]
layer1 = PatternedLayerDefinition(layerDivisions, layer1depth, layerPattern1)

# Reflection and transmission layers
bottomLayer = SemiInfiniteLayerDefinition("bottom material")
topLayer = SemiInfiniteLayerDefinition("top material")

layerStack = [bottomLayer, layer1, layer2, topLayer]


######## DEFINE HARMONICSTRUNCATION ################################################################
# Using 3x3-order only:
M,N = 1,1
harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)

######### DEFINE ANALYSIS ##########################################################################
analysisDefinition = RelativeReflectedTransmittedOrdersAnalysisDefinition()


######## DEFINE SIMULATIONDEFINITION ###############################################################
simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition)


######### RUN SIM #################################################################################
simResults = runSimulation(simulationDefinition)
totalTransmittance = sum(real(simResults[:outputTopRelativeFlux]))
totalReflectance = sum(real(simResults[:outputBottomRelativeFlux]))

# @test isapprox(simResults[:outputBottomRelativeFlux], Rbenchmark, rtol=1e-2)
# @test isapprox(totalReflectance, 0.088768, rtol=1e-2)
# @test isapprox(simResults[:outputTopRelativeFlux], Tbenchmark, rtol=1e-2)
# @test isapprox(totalTransmittance, 0.91123, rtol=1e-2)

######## END OF SHORTEST SIMULATION ############################################################




######## CALCULATE DERIVEDPARAMETERS ###############################################################
derivedParameters = DerivedParameters(simulationDefinition)

# Harmonic orders
@test derivedParameters.harmonicsSet.mnᵢ == [_2VectorInt(1, 1), _2VectorInt(0, 1), _2VectorInt(-1, 1), _2VectorInt(1,0), _2VectorInt(0, 0), _2VectorInt(-1,0), _2VectorInt(1, -1), _2VectorInt(0,-1), _2VectorInt(-1, -1),]
@test derivedParameters.harmonicsSet.indᵢ_mn == Dict(_2VectorInt(1, 1)=>1, _2VectorInt(0,1)=>2, _2VectorInt(-1, 1)=>3, _2VectorInt(1,0)=>4, _2VectorInt(0, 0)=>5, _2VectorInt(-1,0)=>6, _2VectorInt(1, -1)=>7, _2VectorInt(0,-1)=>8, _2VectorInt(-1, -1)=>9)

# k-vector of boundary conditions
@test isapprox(derivedParameters.boundaryConditions.kXY₀/getk₀(simulationDefinition), kBenchmark[X:Y], rtol=1e-3)

# 0-order kXY-vector
#old:
# @test isapprox(derivedParameters.kVectorSet.kᵢ[5]/getk₀(simulationDefinition), kBenchmark[X:Y], rtol=1e-3)
# KNORM
@test isapprox(derivedParameters.kVectorSet.kᵢNorm[5], kBenchmark[X:Y], rtol=1e-3)

# kz components at top and bottom layers
#old:
# kzᵦ = Diagonal( derivedParameters.kzBottom )
# @test isapprox(kzᵦ/getk₀(simulationDefinition), -1*kzᵦbenchmark, rtol=1e-3)
# kzₜ = Diagonal( derivedParameters.kzTop )
# @test isapprox(kzₜ/getk₀(simulationDefinition), kzₜbenchmark, rtol=1e-3)
# Norm
kzᵦNorm = Diagonal( derivedParameters.kzNormBottom )
@test isapprox(kzᵦNorm, kzᵦbenchmark, rtol=1e-3)
# @test isapprox(kzᵦNorm, -1*kzᵦbenchmark, rtol=1e-3)
kzₜNorm = Diagonal( derivedParameters.kzNormTop )
@test isapprox(kzₜNorm, kzₜbenchmark, rtol=1e-3)

# Free space parameters
@test isapprox(derivedParameters.freeSpaceParameters.KzNorm, KzNormBenchmark, rtol=1e-3)
@test isapprox(derivedParameters.freeSpaceParameters.Q, Qbenchmark, rtol=1e-3)
@test isapprox(derivedParameters.freeSpaceParameters.W₀.matrix, Array{ComplexF64,2}(I,(9*2,9*2)), rtol=1e-3)
@test isapprox(derivedParameters.freeSpaceParameters.Λ, λeigenvaluesBenchmark, rtol=1e-3)
@test isapprox(derivedParameters.freeSpaceParameters.V₀.matrix, V₀benchmark, rtol=1e-3)
@test isapprox(derivedParameters.freeSpaceParameters.V₀.matrix, V₀benchmark, rtol=1e-3)
@test isapprox(derivedParameters.boundaryConditions.kXY₀/getk₀(boundaryDefinition.wavenumber), kBenchmark[X:Y], rtol=1e-3)

# K-vector set
# old
# @test isapprox(derivedParameters.kVectorSet.kᵢ[5]/getk₀(derivedParameters.kVectorSet.wavenumber), kBenchmark[X:Y], rtol=1e-3)
# KNORM
@test isapprox(derivedParameters.kVectorSet.kᵢNorm[5], kBenchmark[X:Y], rtol=1e-3)

# Free space parameters.
@test isapprox(derivedParameters.freeSpaceParameters.KzNorm, KzNormBenchmark, rtol=1e-3)
@test isapprox(derivedParameters.freeSpaceParameters.Q, Qbenchmark, rtol=1e-3)
@test isapprox(derivedParameters.freeSpaceParameters.W₀.matrix, Array{ComplexF64,2}(I,(9*2,9*2)), rtol=1e-3)
@test isapprox(derivedParameters.freeSpaceParameters.Λ, λeigenvaluesBenchmark, rtol=1e-3)
@test isapprox(derivedParameters.freeSpaceParameters.V₀.matrix, V₀benchmark, rtol=1e-3)
@test isapprox(derivedParameters.freeSpaceParameters.V₀.matrix, V₀benchmark, rtol=1e-3)


######## CALC INPUT FIELDS #########################################################################

# inputFields = calcInputFields(simulationDefinition, derivedParameters)
inputFields = calcInputFields(derivedParameters.boundaryConditions, derivedParameters.harmonicsSet, derivedParameters.kVectorSet, simulationDefinition.layerStack, simulationDefinition.materialCollection, getWavenumber(simulationDefinition) )
# @test isapprox(inputFields.bottom.modeFields, sourceFields1benchmark, rtol=1e-3)
# @test isapprox(inputFields.top.modeFields, sourceFields2benchmark, rtol=1e-3)


######## CALCULATE GLOBAL SCATTERING MATRIX ########################################################

Sglobal = calcGlobalScatteringMatrix(simulationDefinition, derivedParameters)

_1, _2 = getQuadrantSlices(numHarmonics(derivedParameters.kVectorSet))
# @test isapprox(Sglobal.matrix[_1,_1], SG₁₁benchmark,rtol=1e-3)
# @test isapprox(Sglobal.matrix[_1,_2], SG₁₂benchmark,rtol=1e-3)
# @test isapprox(Sglobal.matrix[_2,_1], SG₂₁benchmark,rtol=1e-3)
# @test isapprox(Sglobal.matrix[_2,_2], SG₂₂benchmark,rtol=1e-3)

######## CALCULATE OUTPUT FIELDS ###############################################################

# inputFields = calcInputFields(simulationDefinition, derivedParameters)
# @test isapprox(inputFields.bottom.modeFields, sourceFields1benchmark,rtol=1e-3 )

outputFields = propagateFields( Sglobal, inputFields, derivedParameters )
# @test isapprox(outputFields.bottom.modeFields, EᵦxyBenchmark, rtol=1e-2)
# @test isapprox(outputFields.top.modeFields, EₜxyBenchmark, rtol=1e-2)


######## CALCULATE DIFFRACTION EFFICIENCIES #######################################################

# Get absolute flux through top and bottom layers
inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux = calcPowerFluxes(inputFields::InputFields, outputFields::OutputFields, derivedParameters::DerivedParameters)

# Get flux relative to input power
inputBottomRelativeFlux, inputTopRelativeFlux, outputBottomRelativeFlux, outputTopRelativeFlux = calcRelativeFluxes(inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux)

# Sum the power fluxes
totalTransmittance = sum(real(outputTopRelativeFlux))
totalReflectance = sum(real(outputBottomRelativeFlux))

@test isapprox(outputBottomRelativeFlux, Rbenchmark, rtol=1e-2)
@test isapprox(totalReflectance, 0.088768, rtol=1e-2)
@test isapprox(outputTopRelativeFlux, Tbenchmark, rtol=1e-2)
@test isapprox(totalTransmittance, 0.91123, rtol=1e-2)



# Short version:
# inputBottomPowerFlux, inputTopPowerFlux, outputBottomPowerFlux, outputTopPowerFlux = calculateReflectedTransmittedOrders(simulationDefinition)


######## END OF SHORTENED INTEGRATION TEST ########################################################








# CALCULATE RI OF TOP AND BOTTOM LAYERS
nᵦ = convert_ϵ2n(calc_ϵ( getMaterial(matCol,"bottom material"), wavenumber))
nₜ = convert_ϵ2n(calc_ϵ( getMaterial(matCol,"top material"), wavenumber))
@test isapprox( nᵦ, 1.4142, rtol = 1e-3)
@test isapprox( nₜ, 3, rtol = 1e-3)
nᵦ = getn( getBoundaryLayer( layerStack, BOTTOM ), matCol, boundaryDefinition.wavenumber )
nₜ = getn( getBoundaryLayer( layerStack, TOP ), matCol, boundaryDefinition.wavenumber )
@test isapprox( nᵦ, 1.4142, rtol = 1e-3)
@test isapprox( nₜ, 3, rtol = 1e-3)

@test isapprox( getk₀(wavenumber), 314.1593, rtol = 1e-3)



# Calculate harmonics set
harmonicsSet = calcHarmonicsSet(harmonicsTruncation)
@test harmonicsSet.mnᵢ == [_2VectorInt(1, 1), _2VectorInt(0, 1), _2VectorInt(-1, 1), _2VectorInt(1,0), _2VectorInt(0, 0), _2VectorInt(-1,0), _2VectorInt(1, -1), _2VectorInt(0,-1), _2VectorInt(-1, -1),]
@test harmonicsSet.indᵢ_mn == Dict(_2VectorInt(1, 1)=>1, _2VectorInt(0,1)=>2, _2VectorInt(-1, 1)=>3, _2VectorInt(1,0)=>4, _2VectorInt(0, 0)=>5, _2VectorInt(-1,0)=>6, _2VectorInt(1, -1)=>7, _2VectorInt(0,-1)=>8, _2VectorInt(-1, -1)=>9)

boundaryConditions = InputByOrderBoundaryConditions( boundaryDefinition, matCol, layerStack )
@test isapprox(boundaryConditions.kXY₀/getk₀(boundaryDefinition.wavenumber), kBenchmark[X:Y], rtol=1e-3)

# Backcalculate kz.  Not used anywhere else.
kzPositive = true
nk₀ = getk₀(wavenumber)*nᵦ
kz² = nk₀^2 - boundaryConditions.kXY₀[X]^2 - boundaryConditions.kXY₀[Y]^2
# kz² = nk₀^2 - k[X]^2 - k[Y]^2
kzIncident = sqrt( kz² )*bool2posNeg(kzPositive)
@test isapprox(kzIncident, 0.70711*getk₀(wavenumber), rtol = 1e-3)
kzIncident = getkz(boundaryConditions.kXY₀[X:Y], nᵦ, wavenumber, kzPositive )
@test isapprox(kzIncident, 0.70711*getk₀(wavenumber), rtol = 1e-3)

kXYZ = kXYtokXYZ(boundaryConditions.kXY₀[X:Y], nᵦ, wavenumber, kzPositive )
@test isapprox(kXYZ/getk₀(wavenumber), kBenchmark, rtol=1e-3)

kzPositive = true
kXYZ = getkXYZ₀(boundaryConditions, layerStack, matCol, kzPositive)
# kXYZ = getkXYZ₀(boundaryDefinition, boundaryConditions, layerStack, matCol)
@test isapprox(kXYZ/getk₀(wavenumber), kBenchmark, rtol=1e-3)





Gvectors = GvectorSet(harmonicsSet,lattice)
kVectorSet = createKVectorSet(boundaryDefinition, boundaryConditions, Gvectors)



# Calculate kx and ky diagonal arrays.  Not used anywhere else.
# old:
# kxArr = Diagonal([kᵢ[X] for kᵢ in kVectorSet.kᵢ])
# kyArr = Diagonal([kᵢ[Y] for kᵢ in kVectorSet.kᵢ])
# @test isapprox(kxArr/getk₀(wavenumber), kxbenchmark, rtol=1e-3)
# @test isapprox(kyArr/getk₀(wavenumber), kybenchmark, rtol=1e-3)
# KNORM
kxArr = Diagonal([kᵢNorm[X] for kᵢNorm in kVectorSet.kᵢNorm])
kyArr = Diagonal([kᵢNorm[Y] for kᵢNorm in kVectorSet.kᵢNorm])
@test isapprox(kxArr, kxbenchmark, rtol=1e-3)
@test isapprox(kyArr, kybenchmark, rtol=1e-3)


# CALCULATE Z-COMPONENTS OF THE K-VECTORS IN THE TOP AND BOTTOM LAYERS
# old:
# kzᵦ = Diagonal( ComplexF64[ conj(sqrt((getk₀(kVectorSet)*nᵦ)^2 - kᵢ[X]^2 - kᵢ[Y]^2))  for kᵢ in kVectorSet.kᵢ] )
# @test isapprox(kzᵦ/getk₀(kVectorSet), -1*kzᵦbenchmark, rtol=1e-3)
# kzₜ = Diagonal( ComplexF64[ conj(sqrt((getk₀(kVectorSet)*nₜ)^2 - kᵢ[X]^2 - kᵢ[Y]^2))  for kᵢ in kVectorSet.kᵢ] )
# @test isapprox(kzₜ/getk₀(kVectorSet), kzₜbenchmark, rtol=1e-3)
# KNORM
# kzᵦ = Diagonal( ComplexF64[ conj(sqrt( nᵦ^2 - kᵢ[X]^2 - kᵢ[Y]^2))  for kᵢ in kVectorSet.kᵢNorm] )
# @test isapprox(kzᵦ, kzᵦbenchmark, rtol=1e-3)
# kzₜ = Diagonal( ComplexF64[ conj(sqrt( nₜ^2 - kᵢ[X]^2 - kᵢ[Y]^2))  for kᵢ in kVectorSet.kᵢNorm] )
# @test isapprox(kzₜ, kzₜbenchmark, rtol=1e-3)

# old
# kzᵦ = Diagonal( calckz(kVectorSet, bottomLayer, matCol, wavenumber) )
# @test isapprox(kzᵦ/getk₀(kVectorSet), -1*kzᵦbenchmark, rtol=1e-3)
# kzₜ = Diagonal( calckz(kVectorSet, topLayer, matCol, wavenumber) )
# @test isapprox(kzₜ/getk₀(kVectorSet), kzₜbenchmark, rtol=1e-3)
# KNORM
kzᵦ = Diagonal( calckzBottom(kVectorSet, bottomLayer, matCol, wavenumber) )
@test isapprox(kzᵦ, kzᵦbenchmark, rtol=1e-3)
kzₜ = Diagonal( calckzTop(kVectorSet, topLayer, matCol, wavenumber) )
@test isapprox(kzₜ, kzₜbenchmark, rtol=1e-3)


freeSpaceParameters = FreeSpaceParameters(derivedParameters.kVectorSet)
@test isapprox(freeSpaceParameters.KzNorm, KzNormBenchmark, rtol=1e-3)
@test isapprox(freeSpaceParameters.Q, Qbenchmark, rtol=1e-3)
@test isapprox(freeSpaceParameters.W₀.matrix, Array{ComplexF64,2}(I,(9*2,9*2)), rtol=1e-3)
@test isapprox(freeSpaceParameters.Λ, λeigenvaluesBenchmark, rtol=1e-3)
@test isapprox(freeSpaceParameters.V₀.matrix, V₀benchmark, rtol=1e-3)
@test isapprox(freeSpaceParameters.V₀.matrix, V₀benchmark, rtol=1e-3)


inputFields = calcInputFields(boundaryConditions, harmonicsSet, kVectorSet, layerStack, matCol, wavenumber)
# inputFields = calcInputFields(boundaryConditions, harmonicsSet, kVectorSet, bottomLayer, topLayer, matCol, wavenumber)
@test isapprox(inputFields.bottom.modeFields, sourceFields1benchmark, rtol=1e-3)
@test isapprox(inputFields.top.modeFields, sourceFields2benchmark, rtol=1e-3)


##### CALCULATE CONVOLUTION MATRICES OF INTERMEDIATE LAYERS
Cϵᵢⱼ1, Cμᵢⱼ1 = calcConvolutionMatrices( layer1, lattice, Gvectors, matCol, wavenumber )
@test Cμᵢⱼ1 ≈ Array{ComplexF64,2}(I,(9,9))
@test isapprox(Cϵᵢⱼ1, Cϵᵢⱼ1benchmark, rtol=1e-1)

Cϵᵢⱼ2, Cμᵢⱼ2 = calcConvolutionMatrices( layer2, lattice, Gvectors, matCol, wavenumber )
@test Cϵᵢⱼ2[1,1] ≈ Complex(6,0)
@test Cμᵢⱼ2[1,1] ≈ Complex(1,0)
@test Cϵᵢⱼ2 ≈ 6*Array{ComplexF64,2}(I,(9,9))
@test Cμᵢⱼ2 ≈ Array{ComplexF64,2}(I,(9,9))

#Calculate inverse convolution matrices for each layer
Cϵᵢⱼ⁻¹1 = inv(Cϵᵢⱼ1)
Cμᵢⱼ⁻¹1 = inv(Cμᵢⱼ1)
Cϵᵢⱼ⁻¹2 = inv(Cϵᵢⱼ2)
Cμᵢⱼ⁻¹2 = inv(Cμᵢⱼ2)



# Initialize global scattering matrix
# Sglobal = initializeGlobalS(numHarmonics(harmonicsSet))
# Don't need to do this separately
Sglobal = initializeGlobalScatteringMatrix( Float64, numHarmonics(harmonicsSet) )
SglobalBenchmark11 = zeros(ComplexF64,(9*2,9*2))
SglobalBenchmark12 = Array{ComplexF64,2}(I,(9*2,9*2))
SglobalBenchmark21 = Array{ComplexF64,2}(I,(9*2,9*2))
SglobalBenchmark22 = zeros(ComplexF64,(9*2,9*2))
SglobalBenchmark = vcat( hcat(SglobalBenchmark11, SglobalBenchmark12),
                        hcat(SglobalBenchmark21, SglobalBenchmark22) )
@test Sglobal.matrix ≈ SglobalBenchmark

# STEP 7: MAIN LOOP
PrecisionType = Float32
prealloc = ScatteringMatrixAllocations{PrecisionType}(numHarmonics(kVectorSet), kVectorSet)

# Layer 1: Patterned
P₁ = calcPmatrixPatterned(prealloc, kVectorSet, Cϵᵢⱼ1, Cϵᵢⱼ⁻¹1, Cμᵢⱼ1, Cμᵢⱼ⁻¹1)
@test isapprox(P₁, P₁benchmark, rtol=1e-3)

Q₁ = calcQmatrixPatterned(prealloc, kVectorSet, Cϵᵢⱼ1, Cϵᵢⱼ⁻¹1, Cμᵢⱼ1, Cμᵢⱼ⁻¹1)
@test isapprox(Q₁, Q₁benchmark, rtol=1e-3)

Ω²₁ = calcΩ²(prealloc, P₁,Q₁)
@test isapprox(Ω²₁, Ω²₁benchmark, rtol=1e-3)


# NOTE: The eigenvalue decomposition is not unique, so values relying on it will not be identical to benchmark, before calculation of the full scattering matrix.
W₁, λ₁ = calcWᵢλᵢ(prealloc, Ω²₁)
@test Ω²₁ * W₁.matrix ≈ W₁.matrix * (λ₁.^2)  # Confirming that it is an eigendecomposition
@test isapprox(Ω²₁benchmark * W₁benchmark, W₁benchmark * (λ₁benchmark.^2), rtol=1e-3)
@test isapprox(Ω²₁benchmark * W₁benchmark, W₁benchmark * (λ₁benchmark^2), rtol=1e-3)

V₁ = calcMagneticEigenvectorsFromQWλ(prealloc, Q₁, W₁, λ₁)

A₁, B₁ = calcAB(prealloc, W₁, derivedParameters.freeSpaceParameters.W₀,V₁, derivedParameters.freeSpaceParameters.V₀)
X₁ = calcX(prealloc, λ₁, kVectorSet.wavenumber, layer1.thickness)

_1, _2 = getQuadrantSlices(numHarmonics(kVectorSet))

S₁ = calcScatteringMatrix_ABX(prealloc, A₁,B₁,X₁)
@test isapprox(S₁.matrix[_1,_1], S1₁₁benchmark, rtol=1e-3)
@test isapprox(S₁.matrix[_1,_2], S1₁₂benchmark, rtol=1e-3)
@test isapprox(S₁.matrix[_2,_1], S1₂₁benchmark, rtol=1e-3)
@test isapprox(S₁.matrix[_2,_2], S1₂₂benchmark, rtol=1e-3)
# sugary:
S₁ = calcScatteringMatrix(prealloc, layer1, matCol, kVectorSet, Gvectors, lattice )
@test isapprox(S₁.matrix[_1,_1], S1₁₁benchmark, rtol=1e-3)
@test isapprox(S₁.matrix[_1,_2], S1₁₂benchmark, rtol=1e-3)
@test isapprox(S₁.matrix[_2,_1], S1₂₁benchmark, rtol=1e-3)
@test isapprox(S₁.matrix[_2,_2], S1₂₂benchmark, rtol=1e-3)



# Layer 2: Unpatterned
mat₂str = layer2.backgroundMaterialName
mat₂ = getMaterial(matCol,mat₂str)
ϵ₂, μ₂ = calc_ϵμ(mat₂,wavenumber)
P₂ = calcPmatrixUnpatterned(prealloc, kVectorSet, ϵ₂, μ₂ )
@test isapprox(P₂, P₂benchmark, rtol=1e-3)
Q₂ = calcQmatrixUnpatterned(prealloc, P₂, ϵ₂, μ₂)
@test isapprox(Q₂, Q₂benchmark, rtol=1e-3)
# Sugary version.  Method specific to uniform layer
P₂, Q₂ = calcPQmatrix(prealloc, layer2, kVectorSet, matCol)
@test isapprox(P₂, P₂benchmark, rtol=1e-3)
@test isapprox(Q₂, Q₂benchmark, rtol=1e-3)
# calcΩ² same for each layer
Ω²₂ = calcΩ²(prealloc, P₂,Q₂)
@test isapprox(Ω²₂, Ω²₂benchmark, rtol=1e-3)
W₂, λ₂ = calcWᵢλᵢ(prealloc, Ω²₂)
V₂ = calcMagneticEigenvectorsFromQWλ(prealloc, Q₂,W₂,λ₂)
# sugary:
V₂ = calcEigenmodesForUniformLayer(prealloc, kVectorSet, layer2, matCol)

# Common components of scattering matrix
A₂ = calcA(prealloc, W₂, derivedParameters.freeSpaceParameters.W₀, V₂, derivedParameters.freeSpaceParameters.V₀)
B₂ = calcB(prealloc, W₂, derivedParameters.freeSpaceParameters.W₀, V₂, derivedParameters.freeSpaceParameters.V₀)
X₂ = calcX(prealloc, λ₂, kVectorSet.wavenumber, layer2.thickness)

S₂ = calcScatteringMatrix_ABX(prealloc, A₂, B₂, X₂)
@test isapprox(S₂.matrix[_1,_1], S2₁₁benchmark, rtol=1e-3)
@test isapprox(S₂.matrix[_1,_2], S2₁₂benchmark, rtol=1e-3)
@test isapprox(S₂.matrix[_2,_1], S2₂₁benchmark, rtol=1e-3)
@test isapprox(S₂.matrix[_2,_2], S2₂₂benchmark, rtol=1e-3)

S₂ = calcScatteringMatrix(prealloc, layer2, matCol, kVectorSet)
@test isapprox(S₂.matrix[_1,_1], S2₁₁benchmark, rtol=1e-3)
@test isapprox(S₂.matrix[_1,_2], S2₁₂benchmark, rtol=1e-3)
@test isapprox(S₂.matrix[_2,_1], S2₂₁benchmark, rtol=1e-3)
@test isapprox(S₂.matrix[_2,_2], S2₂₂benchmark, rtol=1e-3)



# STEP 8: Reflection side scattering matrix
Pᵦ, Qᵦ = calcPQmatrix(prealloc, bottomLayer, kVectorSet, matCol)
@test isapprox(Qᵦ, Qᵦbenchmark, rtol=1e-3)
Ω²ᵦ = calcΩ²(prealloc, Pᵦ,Qᵦ)

# old
# Λᵦ = Array(vcat( hcat(1im*kzᵦ, zeros(ComplexF64,size(kzᵦ)) ),
#            hcat(zeros(ComplexF64,size(kzᵦ )), 1im*kzᵦ) ) / getk₀(kVectorSet))
# KNORM
Λᵦ = Array(vcat( hcat(-1im*kzᵦ, zeros(ComplexF64,size(kzᵦ)) ),
           hcat(zeros(ComplexF64,size(kzᵦ )), -1im*kzᵦ) ) )   # Lecture 7B
@test isapprox(Λᵦ, Λᵦbenchmark, rtol=1e-3)
Λᵦ = calcΛsemiInfiniteBottom(prealloc, Array(kzᵦ), kVectorSet.wavenumber)
@test isapprox(Λᵦ, Λᵦbenchmark, rtol=1e-3)
Wᵦ = derivedParameters.freeSpaceParameters.W₀
# Wᵦeigenmodes = Wᵦ


# Vᵦ = calcEigenmodesFromQλ(Qᵦ,λᵦ)
Vᵦ = calcMagneticEigenvectorsFromQWλ(prealloc, Qᵦ,Wᵦ,Λᵦ)
@test isapprox(Vᵦ.matrix,Vᵦbenchmark,rtol=1e-3)


Aᵦ = calcA_SemiInfinite(prealloc, Wᵦ, derivedParameters.freeSpaceParameters.W₀, Vᵦ, derivedParameters.freeSpaceParameters.V₀)
@test isapprox(Aᵦ, Aᵦbenchmark, rtol=1e-3)
Bᵦ = calcB_SemiInfinite(prealloc, Wᵦ, derivedParameters.freeSpaceParameters.W₀, Vᵦ, derivedParameters.freeSpaceParameters.V₀)
@test isapprox(Bᵦ, Bᵦbenchmark, rtol=1e-3)

# #sugary
Aᵦ, Bᵦ = calcAB_SemiInfinite(prealloc, Wᵦ, derivedParameters.freeSpaceParameters.W₀, Vᵦ, derivedParameters.freeSpaceParameters.V₀)
@test isapprox(Aᵦ, Aᵦbenchmark, rtol=1e-3)
@test isapprox(Bᵦ, Bᵦbenchmark, rtol=1e-3)

Sᵦ = calcScatteringMatrixBottom_AB(prealloc, Aᵦ,Bᵦ)
@test isapprox(Sᵦ.matrix[_1,_1], SR₁₁benchmark, rtol=1e-3)
@test isapprox(Sᵦ.matrix[_1,_2], SR₁₂benchmark, rtol=1e-3)
@test isapprox(Sᵦ.matrix[_2,_1], SR₂₁benchmark, rtol=1e-3)
@test isapprox(Sᵦ.matrix[_2,_2], SR₂₂benchmark, rtol=1e-3)


Sᵦ = calcScatteringMatrixBottom(prealloc, bottomLayer, matCol, kVectorSet)
@test isapprox(Sᵦ.matrix[_1,_1], SR₁₁benchmark, rtol=1e-3)
@test isapprox(Sᵦ.matrix[_1,_2], SR₁₂benchmark, rtol=1e-3)
@test isapprox(Sᵦ.matrix[_2,_1], SR₂₁benchmark, rtol=1e-3)
@test isapprox(Sᵦ.matrix[_2,_2], SR₂₂benchmark, rtol=1e-3)



# STEP 9: Transmission side scattering matrix
Pₜ, Qₜ = calcPQmatrix(prealloc, topLayer, kVectorSet, matCol)
Ω²ₜ = calcΩ²(prealloc, Pₜ,Qₜ)

Λₜ = calcΛsemiInfiniteTop(prealloc, Array(kzₜ), kVectorSet.wavenumber)
@test isapprox(Λₜ,Λₜbenchmark,rtol=1e-3)
Wₜ = derivedParameters.freeSpaceParameters.W₀
# Wₜeigenmodes = Wₜ

Vₜ = calcMagneticEigenvectorsFromQWλ(prealloc, Qₜ,Wₜ,Λₜ)
@test isapprox(Vₜ.matrix,Vₜbenchmark,rtol=1e-3)


Aₜ = calcA_SemiInfinite(prealloc, Wₜ, derivedParameters.freeSpaceParameters.W₀, Vₜ, derivedParameters.freeSpaceParameters.V₀)
Bₜ = calcB_SemiInfinite(prealloc, Wₜ, derivedParameters.freeSpaceParameters.W₀, Vₜ, derivedParameters.freeSpaceParameters.V₀)

Sₜ = calcScatteringMatrixTop_AB(prealloc, Aₜ,Bₜ)
@test isapprox(Sₜ.matrix[_1,_1], ST₁₁benchmark, rtol=1e-3)
@test isapprox(Sₜ.matrix[_1,_2], ST₁₂benchmark, rtol=1e-3)
@test isapprox(Sₜ.matrix[_2,_1], ST₂₁benchmark, rtol=1e-3)
@test isapprox(Sₜ.matrix[_2,_2], ST₂₂benchmark, rtol=1e-3)


Sₜ = calcScatteringMatrixTop(prealloc, topLayer, matCol, kVectorSet)
@test isapprox(Sₜ.matrix[_1,_1], ST₁₁benchmark, rtol=1e-3)
@test isapprox(Sₜ.matrix[_1,_2], ST₁₂benchmark, rtol=1e-3)
@test isapprox(Sₜ.matrix[_2,_1], ST₂₁benchmark, rtol=1e-3)
@test isapprox(Sₜ.matrix[_2,_2], ST₂₂benchmark, rtol=1e-3)



# STEP 10: Global scattering matrix:
Sdevice = S₁⊗S₂
Sglobal = GlobalScatteringMatrix(Sᵦ⊗Sdevice⊗Sₜ)
@test isapprox(Sglobal.matrix[_1,_1], SG₁₁benchmark,rtol=1e-3)
@test isapprox(Sglobal.matrix[_1,_2], SG₁₂benchmark,rtol=1e-3)
@test isapprox(Sglobal.matrix[_2,_1], SG₂₁benchmark,rtol=1e-3)
@test isapprox(Sglobal.matrix[_2,_2], SG₂₂benchmark,rtol=1e-3)


# Put it in terms of a device stack:
Sglobal = calcGlobalScatteringMatrix(layerStack, matCol, kVectorSet, Gvectors, lattice, PrecisionType)
@test isapprox(Sglobal.matrix[_1,_1], SG₁₁benchmark,rtol=1e-3)
@test isapprox(Sglobal.matrix[_1,_2], SG₁₂benchmark,rtol=1e-3)
@test isapprox(Sglobal.matrix[_2,_1], SG₂₁benchmark,rtol=1e-3)
@test isapprox(Sglobal.matrix[_2,_2], SG₂₂benchmark,rtol=1e-3)



# STEP 11: FIELDS

# Not used anywhere else.
# realkXYZInput = kXYtokXYZ(inputMode, nᵦ)
# @test isapprox(realkXYZInput, k, rtol=1e-3)
# @test isapprox(inputMode.kXY, k[X:Y], rtol=1e-3) #works
# # @test isapprox(getk₀(inputMode.wavenumber), getk₀(input), rtol=1e-3)
# realkXYZInput = kXYtokXYZ(inputMode.kXY, nᵦ, wavenumber, inputMode.kzPositive )
# @test isapprox(realkXYZInput, k, rtol=1e-3)

# Not used anywhere else
nₜ = getn(getBoundaryLayer( layerStack, BOTTOM ), matCol, boundaryDefinition.wavenumber)
kXYZ = kXYtokXYZ(boundaryConditions.kXY₀, nₜ, wavenumber, kzPositive )
ŝ, p̂ = calcŝp̂(kXYZ)
@test isapprox(ŝ, [0.5,-0.86603,0], rtol=1e-3)
@test isapprox(p̂, [-0.43301,-0.25,0.86603]*-1, rtol=1e-3) #Adding a -1 multiplier because benchmark appears to have different p̂


# Source fields:
inputFields = calcInputFields(boundaryConditions, harmonicsSet, kVectorSet, layerStack, matCol, wavenumber)
@test isapprox(inputFields.bottom.modeFields, sourceFields1benchmark,rtol=1e-3 )
inputModeCoeff = inputFields2InputCoefficients(inputFields, Wᵦ)
@test isapprox( inputModeCoeff.bottom.modeCoefficients, sourceModeCoeff1benchmark, rtol=1e-3)

S₁₁ = Sglobal.matrix[_1,_1]
S₁₂ = Sglobal.matrix[_1,_2]
S₂₁ = Sglobal.matrix[_2,_1]
S₂₂ = Sglobal.matrix[_2,_2]

# Compute transmission and reflection mode coefficients ("cref" in benchmark)
bottomModeCoeff = S₁₁*inputModeCoeff.bottom.modeCoefficients
topModeCoeff = S₂₁*inputModeCoeff.bottom.modeCoefficients
@test isapprox(bottomModeCoeff, bottomModeCoeffBenchmark, rtol=1e-2)
@test isapprox(topModeCoeff, topModeCoeffBenchmark, rtol=1e-2)


# simpler:
outputModeCoeff = propagateModeCoeff(Sglobal, inputModeCoeff)
@test isapprox(outputModeCoeff.bottom.modeCoefficients, bottomModeCoeffBenchmark, rtol=1e-2)
@test isapprox(outputModeCoeff.top.modeCoefficients, topModeCoeffBenchmark, rtol=1e-2)


# Not used anywhere else.
# get mode for the [0,0] mode.
# inputMode = getInput

# inputMode = Mode(kXYZ, simulationDefinition.wavenumber, A)
# polarizationVector =

# P = fieldSPtoFieldXYZ(inputMode, nᵦ)
# @test isapprox(P,[0.35355-0.30619im,-0.61237-0.17678im,0+0.61237im], rtol=1e-3 )





# Compute bottomected and transmitted fields ("eref" in benchmark)
Eᵦxy = Wᵦ.matrix*outputModeCoeff.bottom.modeCoefficients
Eₜxy = Wₜ.matrix*outputModeCoeff.top.modeCoefficients
@test isapprox(Eᵦxy,EᵦxyBenchmark, rtol=1e-2)
@test isapprox(Eₜxy,EₜxyBenchmark, rtol=1e-2)
# simpler
# outputFields = outputCoefficients2OutputFields(outputModeCoeff, Wᵦ, Wₜ)
outputFields = outputCoefficients2OutputFields(outputModeCoeff, Wₜ)
@test isapprox(outputFields.bottom.modeFields,EᵦxyBenchmark, rtol=1e-2)
@test isapprox(outputFields.top.modeFields,EₜxyBenchmark, rtol=1e-2)


Eᵦx = outputFields.bottom.modeFields[1:numHarmonics(kVectorSet)]
Eᵦy = outputFields.bottom.modeFields[ (numHarmonics(kVectorSet)+1):(2*numHarmonics(kVectorSet))]
# old
# Eᵦz = inv(kzᵦ)*(kVectorSet.Kx*Eᵦx + kVectorSet.Ky*Eᵦy)
# KNORM
Eᵦz = -inv(kzᵦ)*(kVectorSet.KxNorm*Eᵦx + kVectorSet.KyNorm*Eᵦy)
Eₜx = outputFields.top.modeFields[1:numHarmonics(kVectorSet)]
Eₜy = outputFields.top.modeFields[ (numHarmonics(kVectorSet)+1):(2*numHarmonics(kVectorSet))]
# old:
# Eₜz = -inv(kzₜ)*(kVectorSet.Kx*Eₜx + kVectorSet.Ky*Eₜy)
# KNORM
Eₜz = -inv(kzₜ)*(kVectorSet.KxNorm*Eₜx + kVectorSet.KyNorm*Eₜy)
@test isapprox(Eᵦx, Eᵦxbenchmark, rtol=1e-2)
@test isapprox(Eᵦy, Eᵦybenchmark, rtol=1e-2)
@test isapprox(Eᵦz, Eᵦzbenchmark, rtol=1e-2)
@test isapprox(Eₜx, Eₜxbenchmark, rtol=1e-2)
@test isapprox(Eₜy, Eₜybenchmark, rtol=1e-2)
@test isapprox(Eₜz, Eₜzbenchmark, rtol=1e-2)

# old
# bottomFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, kVectorSet, derivedParameters.kzBottom)
# topFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.top, kVectorSet, derivedParameters.kzTop)
# bottomFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.bottom, kVectorSet, derivedParameters.kzBottom)
# topFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.top, kVectorSet, derivedParameters.kzTop)
# KNORM
bottomFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, kVectorSet, derivedParameters.kzNormBottom)
topFieldsOutput = convertFieldSetStackToFieldSetXYZ(outputFields.top, kVectorSet, derivedParameters.kzNormTop)
bottomFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.bottom, kVectorSet, derivedParameters.kzNormBottom)
topFieldsInput = convertFieldSetStackToFieldSetXYZ(inputFields.top, kVectorSet, derivedParameters.kzNormTop)
@test isapprox(bottomFieldsOutput.fields[:,X], Eᵦxbenchmark, rtol=1e-2)
@test isapprox(bottomFieldsOutput.fields[:,Y], Eᵦybenchmark, rtol=1e-2)
@test isapprox(bottomFieldsOutput.fields[:,Z], Eᵦzbenchmark, rtol=1e-2)
@test isapprox(topFieldsOutput.fields[:,X], Eₜxbenchmark, rtol=1e-2)
@test isapprox(topFieldsOutput.fields[:,Y], Eₜybenchmark, rtol=1e-2)
@test isapprox(topFieldsOutput.fields[:,Z], Eₜzbenchmark, rtol=1e-2)

# Step 12: Diffraction efficiencies

# old
# kzₜᵢ = calckz(kVectorSet, topLayer, matCol, wavenumber)
# kzᵦᵢ = calckz(kVectorSet, bottomLayer, matCol, wavenumber)
# KNORM
kzₜᵢ = calckzTop(kVectorSet, topLayer, matCol, wavenumber)
kzᵦᵢ = calckzBottom(kVectorSet, bottomLayer, matCol, wavenumber)

inputPowerFlux = calcPowerFlux(bottomFieldsInput, kzᵦᵢ)
transmittedPowerFlux = calcPowerFlux(topFieldsOutput, kzₜᵢ)
reflectedPowerFlux = calcPowerFlux(bottomFieldsOutput, kzᵦᵢ)
totalInputPowerFlux = abs.(sum(real(inputPowerFlux)))

transmittances = real(transmittedPowerFlux) / totalInputPowerFlux
reflectances = real(reflectedPowerFlux) / totalInputPowerFlux
totalTransmittance = abs(sum(real(transmittances)))
totalReflectance = abs(sum(real(reflectances)))
@show totalInputPowerFlux
@test isapprox(abs.(reflectances), Rbenchmark, rtol=1e-2)
@test isapprox(totalReflectance, 0.088768, rtol=1e-2)
@test isapprox(transmittances, Tbenchmark, rtol=1e-2)
@test isapprox(totalTransmittance, 0.91123, rtol=1e-2)

@test true

end;  # End of test set

println("Completed test set.")





end
