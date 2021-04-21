module SweepTests

using Test


println()
println("Start:")
println()

include("../../src/IncludeKrado.jl")

# Defining baseline simulation

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
addMaterial!(matCol, "Si", importYAMLmaterial(raw"MaterialDatabase\data\main\Si\Schinke.yml"))
addMaterial!(matCol, "SiO2", importYAMLmaterial(raw"MaterialDatabase\data\main\SiO2\Gao.yml"))

# Define lattice
U̅ = [100*nm, 100*nm]
lattice = Lattice(U̅)


# Define layers
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformSi = UniformLayerDefinition(1 * μm, "Si")
uniformAir = UniformLayerDefinition(1 * μm, "Air")
semiInfSi = SemiInfiniteLayerDefinition("Si")
layerStackAirSi = [semiInfAir, uniformAir, uniformSi, semiInfSi]


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



# @testset "Spectral sweep" begin





end;



end
