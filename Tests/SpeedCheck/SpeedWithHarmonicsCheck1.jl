# Measure the speed of simple device over a sweep of # of harmonics

module SpeedWithHarmonicsCheck1

module IntegrationTest3
# Using RCWA-Benchmark-Data-3x3 by Raymond Rumpf from empossible.net



println()
println()
println("Start:")

include("../../IncludeKrado.jl")


# Incident wavevector
λ₀ = 500*nm
wavenumber = WavenumberByλ₀(λ₀)

# Benchmark appears to use nonstandard rotation method.  Here, θ is polar angle (rotation around z-axis) and ϕ is zenith angle (rotation around Y-axis).  ϕ rotation occurs first.
θ = 0 * degrees
ϕ = 0 * degrees
Eₚ = 0.70711
Eₛ = -0.70711im
inputAmplitudes = _2VectorComplex(Eₚ, Eₛ)
boundaryDefinition = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)


# Material Collection
matCol = MaterialCollection()
loadMaterial!(matCol, "MaterialDatabase\\data\\main\\SiO2\\Gao.YML", "SiO2")
loadMaterial!(matCol, "MaterialDatabase\\data\\main\\Si\\Aspnes.YML", "Si")
addMaterial!(matCol, "Air", Material(ConstantPermittivity(1)))
n_glass = 1.4 + 0.001im
ϵ_glass = convert_n2ϵ(n_glass)
addMaterial!(matCol, "glass", Material(ConstantPermittivity(ϵ_glass)))

# Lattice
Lx = 200*nm
U̅ = [Lx, 0]
lattice = Lattice(U̅)

# Layers:
layer1depth = 100*nm  # thickness of layer 1
layer1divisions = 500

# rect = Rectangle([0,0], [0.5*Lx,Inf])
# solid = Solid(rect, "glass")
# solid = Solid(rect, "Si")
# layerPattern = LayerPattern( solid, "Air" )

substrate = SemiInfiniteLayerDefinition("Air")
# patternedLayer = PatternedLayerDefinition(layer1divisions, layer1depth, layerPattern)
uniformLayer = UniformLayerDefinition(layer1depth, "glass")
superstrate = SemiInfiniteLayerDefinition("Air")

# layerStack = [substrate, patternedLayer, superstrate]
layerStack = [substrate, uniformLayer, superstrate]



# Harmonics truncation
M,N = 1,0
harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)

# Analysis method
analysisDefinition = PerformanceAnalysisDefinition()
# analysisDefinition = TransmittanceReflectanceAnalysisDefinition()

# Define sim
simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition)

# Run sim
# simResults = runSimulation(simulationDefinition)



# Preferred method of sweeping over multiple parameters
function doSimulationByHarmonics(M)
    

    local harmonicsTruncation = HarmonicsTruncationByRectangle(M,0)

    local simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition)

    local data = runSimulation(simulationDefinition)

    return data
end

# Mvalues = [1,3,10,30,100, 300]
Mvalues = [200]

sectionNames = ["derivedParameters", "outputFields", "Bottom layer", "S matrix layer 2", "Multiply layer 2", "Top layer"]
sectionNamesForDisplay = unifyStringLengths( sectionNames )
# sectionNames = ["derivedParameters", "inputFields", "outputFields", "Sglobal"]
# maxLength = maximum( [length(sectionName) for sectionName in sectionNames] )
# sectionNamesForDisplay = [ setLength(sectionName, maxLength) for sectionName in sectionNames]
sweepData = map(doSimulationByHarmonics, Mvalues)
println(sweepData)

# timerOutput = sweepData[:timerOutput]
# println(timerOutput)

# println()
# println("TimerOutput data (ms):")
# for iSection in 1:length(sectionNames)
#     print(sectionNamesForDisplay[iSection],"\t")
#     for step in 1:length(Mvalues)
#         # print( round( 1e-6*TimerOutputs.time(sweepData[step][sectionNames[iSection]] ), digits=3), "\t")
#         print( round( 1e-6*TimerOutputs.time(timerOutput[step][sectionNames[iSection]] ), digits=3), "\t")
#     end
#     print("\r\n")
# end

# println(sweepData[length(sweepData)])




end



end
