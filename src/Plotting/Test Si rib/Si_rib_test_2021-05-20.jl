module Si_rib_test
# Testing the match in performance for Si rib
# 100nm pitch, 30% DC, 180nm tall Si on Fused silica.
include("C:\\Krado\\src\\IncludeKrado.jl")

cd(dirname(@__FILE__))

# userMaterialPath = "C:\\UserMaterials\\"
# userMaterialPath = "favoriteMaterials.jl"
include(favoriteMaterials.jl)

# λ₀ = 0.5636*μm
λ₀ = 0.3*μm
wavenumber = WavenumberByλ₀(λ₀)

numλ₀ = 10
λ₀Start = 0.207 * μm
λ₀Stop = 0.4 * μm
λ₀s = LinRange(λ₀Start, λ₀Stop, numλ₀)
wavenumbers = [ WavenumberByλ₀(λ) for λ in λ₀s]



# Benchmark appears to use nonstandard rotation method.  Here, θ is azimuthal angle (rotation around z-axis) and ϕ is zenith angle (rotation around Y-axis).  ϕ rotation occurs first.
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


importFavoriteMaterial!(matCol::MaterialCollection, materialImporters::Dict{<:Any,<:Any}, names)


# Define lattice
# U̅ = [100*nm, 100*nm]
period = 100*nm
U̅ = [period, 0]
lattice = Lattice(U̅)

DC = 0.7

# Define layers
# From test:
# semiInfAir = SemiInfiniteLayerDefinition("Air")
# uniformAir = UniformLayerDefinition(1 * μm, "Air")
# uniformSiO2 = UniformLayerDefinition(1 * μm, "SiO2")
# semiInfSiO2 = SemiInfiniteLayerDefinition("SiO2")
# layerStack = [semiInfAir, uniformAir, uniformSiO2, semiInfSiO2]

# Define actual Stack
substrate = SemiInfiniteLayerDefinition("Air")
# substrate = SemiInfiniteLayerDefinition("Fused_silica")

numDivisions = 1000
ribShape = Rectangle([0,0],[period*DC,1*μm])
rectSolid = Solid(ribShape, "Si_Endura")
ribLayerPattern = LayerPattern( rectSolid, "Air" )
ribLayer = PatternedLayerDefinition(numDivisions, 180*nm, ribLayerPattern)

# ribLayer = UniformLayerDefinition(180*nm, "Si_Endura")  # PROBLEM APPEARS FOR
# ribLayer = UniformLayerDefinition(180*nm, "Fused_silica")
# ribLayer = UniformLayerDefinition(180*nm, "Si")

superstrate = SemiInfiniteLayerDefinition("Air")
layerStack = [ substrate, ribLayer, superstrate ]

# Define Harmonics.
# M,N = 0,0
M,N = 4,4
harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)

analysisDefinition = ZeroOrderModesAnalysisDefinition(FORWARD)

simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinitionNormal, matCol, analysisDefinition)
data = runSimulation(simulationDefinition)
@show abs(data.Tsp[S])^2
@show abs(data.Tsp[P])^2
@show abs(data.Rsp[S])^2
@show abs(data.Rsp[P])^2


plotCrossSection(simulationDefinition, numDivisions, materialPlottingParameters; scale=μm)

# error()

function SiRibTpTs(isTp, wavenumber)

    θ = 0
    ϕ = 0.01*degrees

    if isTp
        Eₚ = 1
        Eₛ = 0
    else
        Eₚ = 0
        Eₛ = 1
    end
    inputAmplitudes = _2VectorComplex(Eₛ, Eₚ)
    boundaryDefinition = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)

    # Define material collection:
    matCol = MaterialCollection()
    addMaterial!(matCol, "Air", Material(ConstantPermittivity(1)))
    addMaterial!(matCol, "Si_Endura", importλnkTextMaterial(userMaterialPath*"Si_Endura.txt"; scale=μm, skipRows=1, delimiter=" "))
    addMaterial!(matCol, "Fused_silica", importλnkTextMaterial(userMaterialPath*"Fused_silica_nonzero_k.txt"; scale=nm, skipRows=3, delimiter="\t"))
    addMaterial!(matCol, "Si", importYAMLmaterial(raw"C:\refractiveindex.info-database\database\data\main\Si\Aspnes.yml"))

    # addMaterial!(matCol, "Si", importYAMLmaterial(raw"..\..\..\MaterialDatabase\data\main\Si\Schinke.yml"))
    # addMaterial!(matCol, "SiO2", importYAMLmaterial(raw"..\..\..\MaterialDatabase\data\main\SiO2\Gao.yml"))


    # Define lattice
    # U̅ = [100*nm, 100*nm]
    period = 100*nm
    DC = 0.3
    # lattice = Lattice(period)
    lattice = Lattice([period,0],[0,period])

    # Define Stack
    substrate = SemiInfiniteLayerDefinition("Air")
    # substrate = SemiInfiniteLayerDefinition("Fused_silica")

    numDivisions = 100
    ribShape = Rectangle([0,0],[period*DC,1*μm])
    rectSolid = Solid(ribShape, "Si_Endura")
    ribLayerPattern = LayerPattern( rectSolid, "Air" )
    ribLayer = PatternedLayerDefinition(numDivisions, 180*nm, ribLayerPattern)

    # ribLayer = UniformLayerDefinition(180*nm, "Si_Endura")  # PROBLEM APPEARS FOR
    # ribLayer = UniformLayerDefinition(180*nm, "Fused_silica")
    # ribLayer = UniformLayerDefinition(180*nm, "Si")

    superstrate = SemiInfiniteLayerDefinition("Air")
    layerStack = [ substrate, ribLayer, superstrate ]

    # Define Harmonics.
    M,N = 0,0
    harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)

    analysisDefinition = ZeroOrderModesAnalysisDefinition(FORWARD)

    simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition)

    results = runSimulation(simulationDefinition)
    return results
end

results = [ [ SiRibTpTs(isTp, wavenumber) for isTp in [true, false]] for wavenumber in wavenumbers]

resultsFilename = "Si_rib_test_2021-04-05.txt"


open(resultsFilename,"w") do io

    println(io, "Wavelength (nm)\tTp\tTs\tRp\tRs")
    for iλ in UnitRange(1,length(wavenumbers))
        resultForλ = results[iλ]
        line = string(getλ₀(wavenumbers[iλ])) * "\t"
        # line *= string(resultForλ[1].totalTransmittance) * "\t" # transmittance reflectance
        # line *= string(resultForλ[2].totalTransmittance) * "\t"
        # line *= string(resultForλ[1].totalReflectance) * "\t"
        # line *= string(resultForλ[2].totalReflectance) * "\t"

        # line *= string(resultForλ[1].outputTopRelativeFlux) * "\t"  # Relative reflectance transmittance orders
        # line *= string(resultForλ[2].outputTopRelativeFlux) * "\t"
        # line *= string(resultForλ[1].outputBottomRelativeFlux) * "\t"
        # line *= string(resultForλ[2].outputBottomRelativeFlux) * "\t"

        line *= string(abs(resultForλ[1].Tsp[P])^2) * "\t"  # zero order modes
        line *= string(abs(resultForλ[2].Tsp[S])^2) * "\t"
        line *= string(abs(resultForλ[1].Rsp[P])^2) * "\t"
        line *= string(abs(resultForλ[2].Rsp[S])^2) * "\t"

        println(io,line)
    end

end
# @show results[1]
# @show results

#Plotting

materialPlottingParameters = Dict{String,PlottingParameters}([
    (
        "Air",
        PlottingParameters(;
            color = [255, 255, 240]./255,
            alpha = 0,
            shade = false,
            lineWidth = 0,
            lineStyle = "None",
        ),
    ),
    (
        "Si_Endura",
        PlottingParameters(;
            color = [75, 100, 50]./255,
            alpha = 0,
            shade = true,
            lineStyle = "-",
            lineColor = "r",
        ),
    ),
    (
        "Fused_silica",
        PlottingParameters(;
            color = [100, 255, 255]./255,
            alpha = 0,
            shade = true,
            lineStyle = "-",
            lineColor = [0.2, 0.7, 0.2],
        ),
    ),
])


# plotHarmonicsSet(simulationDefinition)

# plot Cross-section
# UVstart = [-0.5, 0]
# UVstop = [0.5, 0]
# numDivisions = 50
# XYstart = convertUVtoXY(lattice, UVstart)
# XYstop = convertUVtoXY(lattice, UVstop)
# plotCrossSection(simulationDefinition, XYstart, XYstop, numDivisions, materialPlottingParameters; scale=nm)

# numDivisionsXY = 20
# numDivisionsZ = 30
# positionLineXY = PositionGridXYbyMidpoint( XYstart, XYstop, numDivisionsXY)
# plotCrossSectionNbyArray(simulationDefinition, positionLineXY, numDivisionsZ; scale = μm)

end; #module
