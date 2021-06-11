module FresnelCoefficientsTest
using Test


# println("Fresnel Coefficients - Si and SiO2")
@testset "Fresnel Coefficients - Si and SiO2" begin

include("../../src/IncludeKrado.jl")

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
layerStackSiO2Thin = LayerStack([semiInfAir, uniformAirThin, uniformSiO2Thin, semiInfSiO2])
layerStackSiO2Thick = LayerStack([semiInfAir, uniformAirThick, uniformSiO2Thick, semiInfSiO2])


# Define Harmonics.
M,N = 0,0
harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)

analysisDefinition = ZeroOrderModesAnalysisDefinition(FORWARD)

@testset "SiO2 - visible - uniform - thin - angled" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thin, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSiO2_564nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, rₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, rₛBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[P])^2, tₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[S])^2, tₛBenchmark^2, rtol=1e-3)
end;

@testset "SiO2 - visible - uniform - thick - angled" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thick, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSiO2_564nm, 30*degrees)
# @show(data.Rsp)
@test isapprox(abs(data.Rsp[P])^2, rₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, rₛBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[P])^2, tₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[S])^2, tₛBenchmark^2, rtol=1e-3)
end;

@testset "SiO2 - visible - uniform - thin - angled - S" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thin, harmonicsTruncation, boundaryDefinitionAngledJustS, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSiO2_564nm, 30*degrees)
@test isapprox(abs(data.Rsp[S])^2, rₛBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[S])^2, tₛBenchmark^2, rtol=1e-3)
end;

@testset "SiO2 - visible - uniform - thin - angled - P" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thin, harmonicsTruncation, boundaryDefinitionAngledJustP, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSiO2_564nm, 30*degrees)
@test isapprox(abs(data.Tsp[P])^2, tₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[P])^2, rₚBenchmark^2, rtol=1e-3)
end;

@testset "SiO2 - visible - uniform - thin - normal" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thin, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSiO2_564nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, rₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, rₛBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[P])^2, tₚBenchmark^2, rtol=1e-3)
@test isapprox(abs(data.Tsp[S])^2, tₛBenchmark^2, rtol=1e-3)
end;


@testset "SiO2 - visible - patterned - normal" begin
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirPatt = UniformLayerDefinition(1 * nm, "Air")
uniformSiO2Patt = UniformLayerDefinition(1 * nm, "SiO2_564nm")
semiInfSiO2 = SemiInfiniteLayerDefinition("SiO2_564nm")
layerStackSiO2Patt = LayerStack([semiInfAir, uniformAirPatt, uniformSiO2Patt, semiInfSiO2])

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
θ = 0*degrees  # Doesn't appear to make a difference
ϕ = 30*degrees
boundaryDefinitionAngled = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)
θ = 0*degrees
ϕ = 0*degrees
boundaryDefinitionNormal = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)


# Define layers for Si.  Checking for absorptive material.
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAir = UniformLayerDefinition(1 * μm, "Air")
uniformSi = UniformLayerDefinition(1 * μm, "Si564nm")
semiInfSi = SemiInfiniteLayerDefinition("Si564nm")
layerStackSi = LayerStack([semiInfAir, uniformAir, uniformSi, semiInfSi])


@testset "Si uniform - visible" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi564nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)

simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi564nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;





# Same but the film is a patterned layer
 # This leads to an error when using more than M,N=0,0.  But this should only be a problem for "uniform" films using a patterned geometry
numDivisions = [20,20]
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirPatt = PatternedLayerDefinition(numDivisions, 1 * μm, LayerPattern("Air"))
uniformSiPatt = PatternedLayerDefinition(numDivisions, 1 * μm, LayerPattern("Si564nm"))
semiInfSi = SemiInfiniteLayerDefinition("Si564nm")
layerStackSiPatt = LayerStack([semiInfAir, uniformAirPatt, uniformSiPatt, semiInfSi])

@testset "Si patterned - visible" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSiPatt, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi564nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)

simulationDefinition = SimulationDefinition(lattice, layerStackSiPatt, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi564nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end


@testset "Thick patt Si film wth air substrate - visible" begin
#### Same but the substrate is air, and the film is very thick:
semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirPattThick = PatternedLayerDefinition(numDivisions, 1 * mm, LayerPattern("Air"))
uniformSiPattThick = PatternedLayerDefinition(numDivisions, 1 * mm, LayerPattern("Si564nm"))
layerStackSiPattThick = LayerStack([semiInfAir, uniformAirPattThick, uniformSiPattThick, semiInfAir])


simulationDefinition = SimulationDefinition(lattice, layerStackSiPattThick, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi564nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)

simulationDefinition = SimulationDefinition(lattice, layerStackSiPattThick, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi564nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;

@testset "Jones matrix - reflection and transmission" begin
analysisDefinitionJones = JonesMatrixAnalysisDefinition(FORWARD)
simulationDefinition = SimulationDefinition(lattice, layerStackSiO2Thin, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinitionJones)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSiO2_564nm, 30*degrees)
@test isapprox( abs(data.reflection[S,S]), abs(rₛBenchmark), rtol=1e-3, atol=1e-5)
@test isapprox( abs(data.reflection[P,P]), abs(rₚBenchmark), rtol=1e-3, atol=1e-5)
@test isapprox( abs(data.reflection[S,P]), abs(0), rtol=1e-3, atol=1e-5)
@test isapprox( abs(data.reflection[P,S]), abs(0), rtol=1e-3, atol=1e-5)
@test isapprox( abs(data.transmission[S,S]), abs(tₛBenchmark), rtol=1e-3, atol=1e-5)
@test isapprox( abs(data.transmission[P,P]), abs(tₚBenchmark), rtol=1e-3, atol=1e-5)
@test isapprox( abs(data.transmission[S,P]), abs(0), rtol=1e-3, atol=1e-5)
@test isapprox( abs(data.transmission[P,S]), abs(0), rtol=1e-3, atol=1e-5)

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
layerStackSi = LayerStack([semiInfAir, uniformAir, uniformSi, semiInfSi])

@testset "Si uniform - UV - angled incidence" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi302nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;


@testset "Si uniform - UV - normal incidence" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi302nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;

semiInfAir = SemiInfiniteLayerDefinition("Air")
uniformAirThick = UniformLayerDefinition(1 * mm, "Air")
uniformSiThick = UniformLayerDefinition(1 * mm, "Si302nm")
layerStackSiUniformThick = LayerStack([semiInfAir, uniformAirThick, uniformSiThick, semiInfAir])
@testset "Very thick Si uniform - UV - angled incidence" begin
simulationDefinition = SimulationDefinition(lattice, layerStackSiUniformThick, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi302nm, 30*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;


@testset "Very thick Si uniform - UV - normal incidence" begin
# simulationDefinition = SimulationDefinition(lattice, layerStackSiUniformThick, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
simulationDefinition = SimulationDefinition(lattice, layerStackSiUniformThick, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi302nm, 0*degrees)
@test isapprox(abs(data.Rsp[P])^2, abs(rₚBenchmark)^2, rtol=1e-3)
@test isapprox(abs(data.Rsp[S])^2, abs(rₛBenchmark)^2, rtol=1e-3)
end;



@testset "Jones matrix - reflection" begin
analysisDefinitionJones = JonesMatrixAnalysisDefinition(FORWARD)
simulationDefinition = SimulationDefinition(lattice, layerStackSi, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinitionJones)
data = runSimulation(simulationDefinition)
tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi302nm, 30*degrees)
@test isapprox( abs(data.reflection[S,S]), abs(rₛBenchmark), rtol=1e-3, atol=1e-5)
@test isapprox( abs(data.reflection[P,P]), abs(rₚBenchmark), rtol=1e-3, atol=1e-5)
@test isapprox( abs(data.reflection[S,P]), abs(0), rtol=1e-3, atol=1e-5)
@test isapprox( abs(data.reflection[P,S]), abs(0), rtol=1e-3, atol=1e-5)

end;





end;



@testset "WavelengthSweep" begin
# wavelengths = LinRange(start=0.4, stop=0.7, range=5)
nAir = 1

M,N = 0,0
harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)

wavelengths = LinRange(0.4*μm, 0.7*μm, 5)
wavenumbers = [ WavenumberByλ₀(wavelength) for wavelength in wavelengths]

    θ = 0*degrees
    ϕ = 30*degrees
    Eₚ = 1
    Eₛ = 1
    inputAmplitudes = _2VectorComplex(Eₚ, Eₛ)
    boundaryDefinitionAngled = InputByOrderBoundaryDefinition(wavenumbers[1], θ, ϕ, BOTTOM, inputAmplitudes)

matCol = MaterialCollection()
addMaterial!(matCol, "Air", Material(ConstantPermittivity(1)))
addMaterial!(matCol, "Al2O3", importYAMLmaterial("Tests/MaterialImport/MaterialDatabasePiece/data/main/Al2O3/Boidin.yml"))

# Get n,k:
nAl2O3byλ₀ = [ calc_n( getMaterial(matCol,"Al2O3"), wavenumber ) for wavenumber in wavenumbers ]
# tₛBenchmark, tₚBenchmark, rₛBenchmark, rₚBenchmark = calcFresnelCoefficients(nAir, nSi302nm, 30*degrees)
resultsBenchmark = [ calcFresnelCoefficients(nAir, n, ϕ ) for n in nAl2O3byλ₀]

analysisDefinitionJones = JonesMatrixAnalysisDefinition(FORWARD)

U̅ = [100*nm, 100*nm]
lattice = Lattice(U̅)

semiInfAir = SemiInfiniteLayerDefinition("Air")
# uniformAirThin = UniformLayerDefinition(1 * nm, "Air")
# uniformSiO2Thin = UniformLayerDefinition(1 * nm, "SiO2_564nm")
# uniformAir = UniformLayerDefinition(10 * μm, "Air")
# uniformSiO2Thick = UniformLayerDefinition(10 * μm, "SiO2_564nm")
semiInfAl2O3 = SemiInfiniteLayerDefinition("Al2O3")
# layerStackSiO2Thin = LayerStack([semiInfAir, uniformAirThin, uniformSiO2Thin, semiInfSiO2])
layerStack = LayerStack([semiInfAir, semiInfAl2O3])

simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinitionAngled, matCol, analysisDefinitionJones)

wavenumberSweep = WavenumberSweep(wavenumbers, simulationDefinition)

results = runSweep(wavenumberSweep)

for i_wavenumber in 1:length(wavenumbers)
    @test isapprox( results[i_wavenumber].transmission[S,S], resultsBenchmark[i_wavenumber][1], rtol=1e-3, atol=1e-4 )
    @test isapprox( results[i_wavenumber].transmission[P,P], resultsBenchmark[i_wavenumber][2], rtol=1e-3, atol=1e-4 )
    @test isapprox( results[i_wavenumber].reflection[S,S], resultsBenchmark[i_wavenumber][3], rtol=1e-3, atol=1e-4 )
    @test isapprox( results[i_wavenumber].reflection[P,P], resultsBenchmark[i_wavenumber][4], rtol=1e-3, atol=1e-4 )
end

end




end
