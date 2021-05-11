module SandPtransmission

using Test


println()
println("Start:")

@testset "S and P polarization through Si and SiO2" begin

println()
println("Including Krado")

include("../../src/IncludeKrado.jl")

println("Beginning test")

# Incident wavevector
λ₀ = 0.5*μm
wavenumber = WavenumberByλ₀(λ₀)

# Benchmark appears to use nonstandard rotation method.  Here, θ is azimuthal angle (rotation around z-axis) and ϕ is zenith angle (rotation around Y-axis).  ϕ rotation occurs first.
θ = 0*degrees
ϕ = 30*degrees
Eₚ = 1
Eₛ = 1
inputAmplitudes = _2VectorComplex(Eₚ, Eₛ)
boundaryDefinitionAngled = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)

θ = 0*degrees
ϕ = 0*degrees
Eₚ = 1
Eₛ = 1
inputAmplitudes = _2VectorComplex(Eₚ, Eₛ)
boundaryDefinitionNormal = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)


# Define material collection:
matCol = MaterialCollection()
addMaterial!(matCol, "Air", Material(ConstantPermittivity(1)))
addMaterial!(matCol, "SiO2", importYAMLmaterial(raw"c:\refractiveindex.info-database\database\data\main\SiO2\Gao.yml"))
addMaterial!(matCol, "Si", importYAMLmaterial(raw"C:\refractiveindex.info-database\database\data\main\Si\Aspnes.yml"))

nSi564nm = 4.042 + 0.032im  # doing a conjugate makes extra tests fail
nSi302nm = 5.020 + 3.979im

addMaterial!(matCol, "Si564nm", Material(ConstantPermittivity(convert_n2ϵ(nSi564nm))))
addMaterial!(matCol, "Si302nm", Material(ConstantPermittivity(convert_n2ϵ(nSi302nm))))

# Define lattice
U̅ = [100*nm, 100*nm]
lattice = Lattice(U̅)


# Define layers
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAir = UniformLayerDefinition(1 * μm, "Air")
uniformSiO2 = UniformLayerDefinition(1 * μm, "SiO2")
semiInfSiO2 = SemiInfiniteLayerDefinition("SiO2")
layerStack = [semiInfAir, uniformAir, uniformSiO2, semiInfSiO2]


# Define Harmonics.  # TODO: SOLVE PROBLEM WITH M,N = 4,4 or higher
M,N = 0,0
# M,N = 4,4
harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)

analysisDefinition = ZeroOrderModesAnalysisDefinition(FORWARD)

@testset "SiO2 - visible" begin
simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmark = 0.023607
RsBenchmark = 0.054686
@test isapprox(abs(data.Rsp[P])^2, RpBenchmark, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmark, rtol=1e-3) # Both S and P for normal incidence


simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmarkNormal = 0.037664
RsBenchmarkNormal = 0.037664
@test isapprox(abs(data.Rsp[P])^2, RpBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
end;


### SI:


λ₀ = 0.5636*μm
wavenumber = WavenumberByλ₀(λ₀)
θ = 0.001*degrees
# θ = 0*degrees  # Doesn't appear to make a difference
ϕ = 30*degrees
boundaryDefinitionAngled = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)
θ = 0.001*degrees
ϕ = 0.001*degrees
# θ = 0*degrees
# ϕ = 0*degrees
boundaryDefinitionNormal = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)


# Define layers for Si.  Checking for absorptive material.
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAir = UniformLayerDefinition(1 * μm, "Air")
uniformSi = UniformLayerDefinition(1 * μm, "Si564nm")
semiInfSi = SemiInfiniteLayerDefinition("Si564nm")

layerStackSi = [semiInfAir, uniformAir, uniformSi, semiInfSi]

@testset "Si uniform - visible" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmark = 0.31168
RsBenchmark = 0.41586
@test isapprox(abs(data.Rsp[P])^2, RpBenchmark, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmark, rtol=1e-3) # Both S and P for normal incidence


simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmarkNormal = 0.36404
RsBenchmarkNormal = 0.36404
@test isapprox(abs(data.Rsp[P])^2, RpBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
end





# Same but the film is a patterned layer
numDivisions = [20,20]
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirPatt = PatternedLayerDefinition(numDivisions, 1 * μm, LayerPattern("Air"))
uniformSiPatt = PatternedLayerDefinition(numDivisions, 1 * μm, LayerPattern("Si564nm"))
semiInfSi = SemiInfiniteLayerDefinition("Si564nm")
layerStackSiPatt = [semiInfAir, uniformAirPatt, uniformSiPatt, semiInfSi]


@testset "Si patterned - visible" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSiPatt, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmark = 0.31168
RsBenchmark = 0.41586
@test isapprox(abs(data.Rsp[P])^2, RpBenchmark, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmark, rtol=1e-3) # Both S and P for normal incidence


simulationDefinition = SimulationDefinition(lattice, layerStackSiPatt, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmarkNormal = 0.36404
RsBenchmarkNormal = 0.36404
@test isapprox(abs(data.Rsp[P])^2, RpBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
end



@testset "Thick patt Si film wth air substrate - visible" begin
#### Same but the substrate is air, and the film is very thick:
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirPattThick = PatternedLayerDefinition(numDivisions, 1 * mm, LayerPattern("Air"))
uniformSiPattThick = PatternedLayerDefinition(numDivisions, 1 * mm, LayerPattern("Si564nm"))
layerStackSiPattThick = [semiInfAir, uniformAirPattThick, uniformSiPattThick, semiInfAir]


simulationDefinition = SimulationDefinition(lattice, layerStackSiPattThick, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmark = 0.31168
RsBenchmark = 0.41586
@test isapprox(abs(data.Rsp[P])^2, RpBenchmark, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmark, rtol=1e-3) # Both S and P for normal incidence


simulationDefinition = SimulationDefinition(lattice, layerStackSiPattThick, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmarkNormal = 0.36404
RsBenchmarkNormal = 0.36404
@test isapprox(abs(data.Rsp[P])^2, RpBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
end;





# Setting wavelength to 300nm:
λ₀ = 0.3024*μm
wavenumber = WavenumberByλ₀(λ₀)
θ = 0.001*degrees
ϕ = 30*degrees
boundaryDefinitionAngled = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)
θ = 0.001*degrees
ϕ = 0.001*degrees
boundaryDefinitionNormal = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)

# Define layers
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAir = UniformLayerDefinition(1 * μm, "Air")
uniformSi = UniformLayerDefinition(1 * μm, "Si302nm")
semiInfSi = SemiInfiniteLayerDefinition("Si302nm")
layerStackSi = [semiInfAir, uniformAir, uniformSi, semiInfSi]

@testset "Si uniform - UV" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmark = 0.31168
RsBenchmark = 0.41586
@test isapprox(abs(data.Rsp[P])^2, RpBenchmark, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmark, rtol=1e-3) # Both S and P for normal incidence


simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmarkNormal = 0.61439
RsBenchmarkNormal = 0.61439
@test isapprox(abs(data.Rsp[P])^2, RpBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
end;



@testset " Very thick Si uniform - UV" begin
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirThick = UniformLayerDefinition(1 * mm, "Air")
uniformSiThick = UniformLayerDefinition(1 * mm, "Si")
layerStackSiUniformThick = [semiInfAir, uniformAirThick, uniformSiThick, semiInfAir]


simulationDefinition = SimulationDefinition(lattice, layerStackSiUniformThick, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmark = 0.31168
RsBenchmark = 0.41586
@test isapprox(abs(data.Rsp[P])^2, RpBenchmark, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmark, rtol=1e-3) # Both S and P for normal incidence


simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmarkNormal = 0.61439
RsBenchmarkNormal = 0.61439
@test isapprox(abs(data.Rsp[P])^2, RpBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
@test isapprox(abs(data.Rsp[S])^2, RsBenchmarkNormal, rtol=1e-3) # Both S and P for normal incidence
end;






end;

end
