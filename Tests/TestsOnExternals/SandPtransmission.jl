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
inputAmplitudesJustS = _2VectorComplex(Eₛ, 0)
inputAmplitudesJustP = _2VectorComplex(0, Eₚ)
boundaryDefinitionAngled = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)

boundaryDefinitionAngledJustS = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudesJustS)
boundaryDefinitionAngledJustP = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudesJustP)

θ = 0*degrees
ϕ = 0*degrees
Eₚ = 1
Eₛ = 1
inputAmplitudes = _2VectorComplex(Eₚ, Eₛ)
boundaryDefinitionNormal = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)

# Define material collection:
matCol = MaterialCollection()
addMaterial!(matCol, "Air", Material(ConstantPermittivity(1)))

nAir = 1
nSiO2_564nm = 1.478503
nSi564nm = 4.042 + 0.032im
nSi302nm = 5.020 + 3.979im

addMaterial!(matCol, "SiO2_564nm", Material(ConstantPermittivity(convert_n2ϵ(nSiO2_564nm))))
addMaterial!(matCol, "Si564nm", Material(ConstantPermittivity(convert_n2ϵ(nSi564nm))))
addMaterial!(matCol, "Si302nm", Material(ConstantPermittivity(convert_n2ϵ(nSi302nm))))

# Define lattice
U̅ = [100*nm, 100*nm]
lattice = Lattice(U̅)


# Define layers
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirThin = UniformLayerDefinition(1 * nm, "Air")
uniformSiO2Thin = UniformLayerDefinition(1 * nm, "SiO2_564nm")
uniformAirThick = UniformLayerDefinition(10 * μm, "Air")
uniformSiO2Thick = UniformLayerDefinition(10 * μm, "SiO2_564nm")
semiInfSiO2 = SemiInfiniteLayerDefinition("SiO2_564nm")
layerStackSiO2Thin = [semiInfAir, uniformAirThin, uniformSiO2Thin, semiInfSiO2]
layerStackSiO2Thick = [semiInfAir, uniformAirThick, uniformSiO2Thick, semiInfSiO2]


# Define Harmonics.  # TODO: SOLVE PROBLEM WITH M,N = 4,4 or higher
M,N = 0,0
# M,N = 4,4
harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)

analysisDefinition = ZeroOrderModesAnalysisDefinition(FORWARD)

@testset "SiO2 - visible" begin
println("SiO2 - visible - angled incidence - thin")
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thin, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSiO2_564nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, rₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, rₛBenchmark^2, rtol=1e-3)


println("SiO2 - visible - angled incidence - thick")
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thick, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSiO2_564nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, rₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, rₛBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[P])^2, tₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[S])^2, tₛBenchmark^2, rtol=1e-3)


println("SiO2 - visible - angled incidence - just S")
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thin, harmonicsTruncation, boundaryDefinitionAngledJustS, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
@test isapprox(abs(data.Rsp[S])^2, rₛBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[S])^2, tₛBenchmark^2, rtol=1e-3)

println("SiO2 - visible - angled incidence - just P")
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thin, harmonicsTruncation, boundaryDefinitionAngledJustP, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
@test isapprox(abs(data.Tsp[P])^2, tₚBenchmark^2, rtol=1e-3)

println("SiO2 - visible - normal incidence")
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thin, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSiO2_564nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, rₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, rₛBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[P])^2, tₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[S])^2, tₛBenchmark^2, rtol=1e-3)
end;


@testset "SiO2 - visible - patterned - normal" begin
println("SiO2 - visible - patterned - normal")
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirPatt = UniformLayerDefinition(1 * nm, "Air")
uniformSiO2Patt = UniformLayerDefinition(1 * nm, "SiO2_564nm")
semiInfSiO2 = SemiInfiniteLayerDefinition("SiO2_564nm")
layerStackSiO2Patt = [semiInfAir, uniformAirPatt, uniformSiO2Patt, semiInfSiO2]

simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Patt, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSiO2_564nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, rₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, rₛBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[P])^2, tₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[S])^2, tₛBenchmark^2, rtol=1e-3)
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
# layerStackSi = [semiInfAir, uniformAir,  semiInfSi]
# layerStackSi = [semiInfAir, uniformSi,  semiInfSi]

@testset "Si uniform - visible" begin
println("Si uniform - visible - angled incidence")
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi564nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)

println("Si uniform - visible - normal incidence")
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi564nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;





# Same but the film is a patterned layer
numDivisions = [20,20]
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirPatt = PatternedLayerDefinition(numDivisions, 1 * μm, LayerPattern("Air"))
uniformSiPatt = PatternedLayerDefinition(numDivisions, 1 * μm, LayerPattern("Si564nm"))
semiInfSi = SemiInfiniteLayerDefinition("Si564nm")

layerStackSiPatt = [semiInfAir, uniformAirPatt, uniformSiPatt, semiInfSi]
# layerStackSiPatt = [semiInfAir, uniformAirPatt, semiInfSi]
# layerStackSiPatt = [semiInfAir, uniformSiPatt, semiInfSi]

@testset "Si patterned - visible" begin
println("Si patterned - visible - angled incidence")
simulationDefinition = SimulationDefinition(lattice, layerStackSiPatt, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi564nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)

println("Si patterned - visible - normal incidence")
simulationDefinition = SimulationDefinition(lattice, layerStackSiPatt, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi564nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end

# error()

@testset "Thick patt Si film wth air substrate - visible" begin
println("Thick patt Si film wth air substrate - visible - angled incidence")
#### Same but the substrate is air, and the film is very thick:
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirPattThick = PatternedLayerDefinition(numDivisions, 1 * mm, LayerPattern("Air"))
uniformSiPattThick = PatternedLayerDefinition(numDivisions, 1 * mm, LayerPattern("Si564nm"))
layerStackSiPattThick = [semiInfAir, uniformAirPattThick, uniformSiPattThick, semiInfAir]


simulationDefinition = SimulationDefinition(lattice, layerStackSiPattThick, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmark = 0.31168
RsBenchmark = 0.41586
@test isapprox(abs(data.Rsp[P])^2, RpBenchmark, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, RsBenchmark, rtol=1e-3)

println("Thick patt Si film wth air substrate - visible - normal incidence")
simulationDefinition = SimulationDefinition(lattice, layerStackSiPattThick, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)

RpBenchmarkNormal = 0.36404
RsBenchmarkNormal = 0.36404
@test isapprox(abs(data.Rsp[P])^2, RpBenchmarkNormal, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, RsBenchmarkNormal, rtol=1e-3)
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

@testset "Si uniform - UV - angled incidence" begin
println("Si uniform - UV - angled incidence")
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi302nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;


@testset "Si uniform - UV - normal incidence" begin
println("Si uniform - UV - normal incidence")
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi302nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;


@testset " Very thick Si uniform - UV - angled incidence" begin
println(" Very thick Si uniform - UV - angled incidence")
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirThick = UniformLayerDefinition(1 * mm, "Air")
uniformSiThick = UniformLayerDefinition(1 * mm, "Si302nm")
layerStackSiUniformThick = [semiInfAir, uniformAirThick, uniformSiThick, semiInfAir]


simulationDefinition = SimulationDefinition(lattice, layerStackSiUniformThick, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi302nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;


@testset " Very thick Si uniform - UV - normal incidence" begin
println(" Very thick Si uniform - UV - normal incidence")
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi302nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;






end;

end
