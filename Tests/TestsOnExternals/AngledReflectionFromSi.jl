module AngledReflectionFromSi

using Test


println()
println("Start:")

@testset "Testing the angled reflection from Si" begin

println()
println("Including Krado")

include("../../src/IncludeKrado.jl")

println("Beginning test")

# Incident wavevector
λ₀ = 0.5*μm
wavenumber = WavenumberByλ₀(λ₀)

# Benchmark appears to use nonstandard rotation method.  Here, θ is azimuthal angle (rotation around z-axis) and ϕ is zenith angle (rotation around Y-axis).  ϕ rotation occurs first.
θ = 0
ϕ = 0
# Only doing Rp
Eₚ = 1
Eₛ = 0
inputAmplitudes = _2VectorComplex(Eₚ, Eₛ)
boundaryDefinitionForward = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)
boundaryDefinitionBackward = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, TOP, inputAmplitudes)

# Define material collection:
matCol = MaterialCollection()
addMaterial!(matCol, "Air", Material(ConstantPermittivity(1)))
# addMaterial!(matCol, "Air", Material(ConstantPermittivity(1.00027897^2)))
# addMaterial!(matCol, "Air", importYAMLmaterial(raw"MaterialDatabase\data\other\mixed gases\air\Ciddor.yml"))
addMaterial!(matCol, "Si", importYAMLmaterial(raw"..\..\..\MaterialDatabase\data\main\Si\Schinke.yml"))
addMaterial!(matCol, "SiO2", importYAMLmaterial(raw"..\..\..\MaterialDatabase\data\main\SiO2\Gao.yml"))
# addMaterial!(matCol, "SiO2", Material(ConstantPermittivity(1.4^2))


# Define lattice
U̅ = [100*nm, 100*nm]
lattice = Lattice(U̅)


# Define layers
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAir = UniformLayerDefinition(1 * μm, "Air")
uniformSi = UniformLayerDefinition(1 * μm, "Si")
semiInfSi = SemiInfiniteLayerDefinition("Si")
layerStackAirAirSiSi = [semiInfAir, uniformAir, uniformSi, semiInfSi]
layerStackAirAirSi = [semiInfAir, uniformAir, semiInfSi]
layerStackAirSiSi = [semiInfAir, uniformSi, semiInfSi]
layerStackAirSi = [semiInfAir, semiInfSi]
layerStackSiSiAirAir = [semiInfSi, uniformSi, uniformAir, semiInfAir]

semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAir = UniformLayerDefinition(1 * μm, "Air")
uniformSiO2 = UniformLayerDefinition(1 * μm, "SiO2")
semiInfSiO2 = SemiInfiniteLayerDefinition("SiO2")
layerStackAirAirSiO2SiO2 = [semiInfAir, uniformAir, uniformSiO2, semiInfSiO2]
layerStackSiO2SiO2AirAir = [semiInfSiO2, uniformSiO2, uniformAir, semiInfAir]


numDivisions = [10,10]
patternedAirLayer = PatternedLayerDefinition(numDivisions, 1 * μm, LayerPattern("Air"))
patternedSiO2Layer = PatternedLayerDefinition(numDivisions, 1 * μm, LayerPattern("SiO2"))
patternedSiLayer = PatternedLayerDefinition(numDivisions, 1 * μm, LayerPattern("Si"))
patternedLayerStackAirAirSiO2SiO2 = [ semiInfAir, patternedAirLayer, patternedSiO2Layer, semiInfSiO2]
# patternedLayerStackAirAirSiO2SiO2 = [ semiInfAir, patternedAirLayer, patternedSiO2Layer, topLayer]
patternedLayerStackSiO2SiO2AirAir = [ semiInfSiO2, patternedSiO2Layer, patternedAirLayer, semiInfAir]
patternedLayerStackAirAirSi = [ semiInfAir, patternedAirLayer, semiInfSi]
patternedLayerStackAirSiSi = [ semiInfAir, patternedSiLayer, semiInfSi]
patternedLayerStackAirAirSiSi = [ semiInfAir, patternedAirLayer, patternedSiLayer, semiInfSi]
# patternedlayerStackSiSiAirAir = [ topLayer, patternedSiLayer, patternedAirLayer, semiInfAir]


# Define Harmonics.
M,N = 1,0
harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)

analysisDefinitionForward = TransmittanceReflectanceAnalysisDefinition(FORWARD)
analysisDefinitionBackward = TransmittanceReflectanceAnalysisDefinition(BACKWARD)

@testset "Normal incidence" begin

# Air to Si. Normal
println()
println("Air, Si.  Bottom incidence.")
simulationDefinition = SimulationDefinition(lattice, layerStackAirSi, harmonicsTruncation, boundaryDefinitionForward, matCol, analysisDefinitionForward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.38676, rtol=1e-3) # Both S and P for normal incidence

# Air to Si. Normal
println()
println("Air, Air, Si.  Bottom incidence.")
simulationDefinition = SimulationDefinition(lattice, layerStackAirAirSi, harmonicsTruncation, boundaryDefinitionForward, matCol, analysisDefinitionForward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.38676, rtol=1e-3)

println()
println("Air, Si, Si.  Bottom incidence.")
simulationDefinition = SimulationDefinition(lattice, layerStackAirSiSi, harmonicsTruncation, boundaryDefinitionForward, matCol, analysisDefinitionForward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.38676, rtol=1e-3) # Both S and P for normal incidence

# Air to SiO2. Normal
println()
println("Air, Air, SiO2, SiO2.  Bottom incidence.")

simulationDefinition = SimulationDefinition(lattice, layerStackAirAirSiO2SiO2, harmonicsTruncation, boundaryDefinitionForward, matCol, analysisDefinitionForward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.037664, rtol=1e-2) # Both S and P for normal incidence

# SiO2 to Air. Injection inside SiO2
println()
println("SiO2, Air.  Bottom incidence.  Injection inside SiO2")

simulationDefinition = SimulationDefinition(lattice, layerStackSiO2SiO2AirAir, harmonicsTruncation, boundaryDefinitionForward, matCol, analysisDefinitionForward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.037664, rtol=1e-2) # Both S and P for normal incidence


# Air to SiO2. Reversed substrate and superstrate
println()
println("SiO2, SiO2, Air, Air.  Boundary backward. Analysis backward.")

simulationDefinition = SimulationDefinition(lattice, layerStackSiO2SiO2AirAir, harmonicsTruncation, boundaryDefinitionBackward, matCol, analysisDefinitionBackward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.037664, rtol=1e-2) # Both S and P for normal incidence


# SiO2 to Air. Reversed substrate and superstrate. Injection inside SiO2
println()
println("SiO2, Air.  Top incidence.")

simulationDefinition = SimulationDefinition(lattice, layerStackAirAirSiO2SiO2, harmonicsTruncation, boundaryDefinitionBackward, matCol, analysisDefinitionBackward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.037664, rtol=1e-2) # Both S and P for normal incidence


# Air to Si. Reversed substrate and superstrate
println()
println("Si, Air.  Top incidence.")

simulationDefinition = SimulationDefinition(lattice, layerStackSiSiAirAir, harmonicsTruncation, boundaryDefinitionBackward, matCol, analysisDefinitionBackward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.38676, rtol=1e-2) # Both S and P for normal incidence



# Air to SiO2. Normal. Using patterned layer.
println()
println("Air, SiO2.  Bottom incidence.  Using patterned layer.")

simulationDefinition = SimulationDefinition(lattice, patternedLayerStackAirAirSiO2SiO2, harmonicsTruncation, boundaryDefinitionForward, matCol, analysisDefinitionForward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.037664, rtol=1e-2) # Both S and P for normal incidence

# SiO2 to Air. Normal. Using patterned layer. SiO2 incidence.
println()
println("Air, SiO2.  Bottom incidence.  Using patterned layer.")

simulationDefinition = SimulationDefinition(lattice, patternedLayerStackAirAirSiO2SiO2, harmonicsTruncation, boundaryDefinitionBackward, matCol, analysisDefinitionBackward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.037664, rtol=1e-2) # Both S and P for normal incidence


# Air to Si. Normal. Using patterned layer.
println()
println("Air, Si.  Bottom incidence.  Using patterned layer.")

simulationDefinition = SimulationDefinition(lattice, patternedLayerStackAirAirSiSi, harmonicsTruncation, boundaryDefinitionForward, matCol, analysisDefinitionForward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.38676, rtol=1e-2) # Both S and P for normal incidence

# Air to Si. Normal. Using patterned layer.
println()
println("Air, Si.  Bottom incidence.  Using patterned layer.")

simulationDefinition = SimulationDefinition(lattice, patternedLayerStackAirAirSi, harmonicsTruncation, boundaryDefinitionForward, matCol, analysisDefinitionForward)
data = runSimulation(simulationDefinition)

@test isapprox(data.totalReflectance, 0.38676, rtol=1e-2) # Both S and P for normal incidence



println()
println("Air, Si.  Bottom incidence.  Using patterned layer. ZeroOrderModes")

zeroOrderAnalysisDefinition = ZeroOrderModesAnalysisDefinition(FORWARD)

simulationDefinition = SimulationDefinition(lattice, patternedLayerStackAirAirSi, harmonicsTruncation, boundaryDefinitionForward, matCol, zeroOrderAnalysisDefinition)
data = runSimulation(simulationDefinition)

totalReflectance = abs(data.Rsp[S])^2 + abs(data.Rsp[P])^2
@test isapprox(totalReflectance, 0.38676, rtol=1e-2) # Both S and P for normal incidence




end;

@testset "Absorptive film" begin



    #
    # simulationDefinition = SimulationDefinition(lattice, patternedLayerStackAirAirSiO2SiO2, harmonicsTruncation, boundaryDefinitionBackward, matCol, analysisDefinitionBackward)
    # data = runSimulation(simulationDefinition)
    #
    # @test isapprox(data.totalReflectance, 0.037664, rtol=1e-3) # Both S and P for normal incidence


end;

end;


end
#### SWEEP OVER WAVELENGTH ####

# Create generation function for Simulation Definition.
# Parameters is an iterable
# function doSimulationByWavelength(wavelength)
#
#     local wavenumber = WavenumberByλ₀(wavelength)
#     local boundaryDefinition = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)
#
#     simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition)
#
#     data = runSimulation(simulationDefinition)
#
#     return data
# end
