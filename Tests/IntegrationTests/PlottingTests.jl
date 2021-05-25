# PLotting tests:

module TestModule

using Test

include("../../src/IncludeKrado.jl")

@testset "Plotting Tests" begin

    println()
    println("Plotting tests:")


    # Incident wavevector
    wavenumber = WavenumberByλ₀(500 * nm)

    Escale = 1

    # Benchmark appears to use nonstandard rotation method.  Here, θ is azimuthal angle (rotation around z-axis) and ϕ is zenith angle (rotation around Y-axis).  ϕ rotation occurs first.
    θ = 0 * degrees
    ϕ = 30 * degrees
    # Eₚ = 0.70711
    # Eₛ = -0.70711
    inputAmplitudesCircularPolarization = _2VectorComplex(√(2), √(2) * 1im)
    inputAmplitudesSpolarization = _2VectorComplex(1, 0)
    inputAmplitudesSpolarizationWeak = _2VectorComplex(1, 0) * 0.5
    inputAmplitudesPpolarization = _2VectorComplex(0, 1)

    mainHarmonicOrder = _2VectorInt(0, 0)
    isTop = false
    Abyϖbottom = Dict{_2VectorInt,_2VectorComplex}()
    Abyϖtop = Dict{_2VectorInt,_2VectorComplex}()

    # Abyϖtop[_2VectorInt(1, 0)] = inputAmplitudesPpolarization  #Evanescent
    # Abyϖtop[_2VectorInt(0, 1)] = inputAmplitudesSpolarization  #Evanescent
    # Abyϖtop[_2VectorInt(-2, -2)] = inputAmplitudesSpolarizationWeak
    Abyϖbottom[_2VectorInt(0, 0)] = inputAmplitudesCircularPolarization

    boundaryDefinition = InputByOrderBoundaryDefinition(
        wavenumber,
        θ,
        ϕ,
        mainHarmonicOrder,
        isTop,
        Abyϖbottom,
        Abyϖtop,
    )


    # Define lattice
    pitch = 1 * μm
    U̅ = [pitch, 0]
    V̅ = [pitch / 2, pitch * sqrt(3) / 2]
    lattice = Lattice(U̅, V̅; originOffsetUV = [-0.5, -0.5])

    # Define stack
    layer1 = UniformLayerDefinition(100 * nm, "Ag")
    layer2solids = [Solid(Rectangle([0, 0], [0.518 * μm, 0.324 * μm]), "Al2O3")]
    layer2 =
        PatternedLayerDefinition([100, 100], 30 * nm, LayerPattern(layer2solids, "Air"))

    layer3solids = [Solid(Circle([0.1 * μm, 0], 0.15 * μm), "Ag")]
    layer3 =
        PatternedLayerDefinition([100, 100], 100 * nm, LayerPattern(layer3solids, "Al2O3"))
    layer4solids = [Solid(
        Polygon( [[0, 0], [0.1, 0], [-0.05, 0.05], [0, -0.1]] * μm, offset=[0.1 * μm, 0]),
        "Ag",
    )]
    layer4 =
        PatternedLayerDefinition([100, 100], 100 * nm, LayerPattern(layer4solids, "Air"))

    layer5circs = [Circle([0.1 * μm, 0], 0.2 * μm), Circle([-0.1 * μm, 0], 0.2 * μm)]
    layer5solids = [Solid(UnionShape(layer5circs), "Ag")]
    layer5 =
        PatternedLayerDefinition([100, 100], 100 * nm, LayerPattern(layer5solids, "Air"))

    layer6solids = [Solid(IntersectionShape(layer5circs), "Ag")]
    layer6 =
        PatternedLayerDefinition([200, 200], 100 * nm, LayerPattern(layer6solids, "Air"))

    layer7shapes = [
        Circle([0.1 * μm, 0], 0.2 * μm),
        Circle([-0.1 * μm, 0], 0.2 * μm),
        Rectangle([0, 0], [0.2 * μm, 0.4 * μm]),
    ]
    layer7solids = [Solid(DifferenceShape(layer7shapes), "Ag")]
    layer7 =
        PatternedLayerDefinition([200, 200], 200 * nm, LayerPattern(layer7solids, "Air"))

    layer8shapes = [Circle([0.1 * μm, 0], 0.2 * μm), Circle([-0.1 * μm, 0], 0.2 * μm)]
    layer8baseShape = Rectangle([0, 0], [0.2 * μm, 0.45 * μm])
    layer8solids = [Solid(SubtractionShape(layer8baseShape, layer8shapes), "Ag")]
    layer8 =
        PatternedLayerDefinition([200, 200], 100 * nm, LayerPattern(layer8solids, "Air"))
    layer9 = UniformLayerDefinition(100 * nm, "Ag")


    substrateLayer = SemiInfiniteLayerDefinition("Al2O3")
    superstrateLayer = SemiInfiniteLayerDefinition("Air")

    layerStack = [
        substrateLayer,
        layer1,
        layer2,
        layer3,
        layer4,
        layer5,
        layer6,
        layer7,
        layer8,
        layer9,
        superstrateLayer,
    ]


    # Import materials:
    matCol = MaterialCollection()
    addMaterial!(matCol,"Air", Material( ConstantPermittivity(1) ) )
    addMaterial!(matCol,"Al2O3", Material( ConstantPermittivity(1.7) ) )
    addMaterial!(matCol,"Ag", Material( ConstantPermittivity(0.5 + 3.9im) ) )


    # Define Harmonics.
    M, N = 3, 3
    # M,N = 10,10
    harmonicsTruncation = HarmonicsTruncationByRectangle(M, N)
    # M,N = 20,10
    # harmonicsTruncation = SuperellipseHarmonicsTruncation(M,N, γSUPERELLIPSE["circle"])

    analysisDefinition = TransmittanceReflectanceAnalysisDefinition(FORWARD)

    simulationDefinition = SimulationDefinition(
        lattice,
        layerStack,
        harmonicsTruncation,
        boundaryDefinition,
        matCol,
        analysisDefinition,
    )


    # Define what materials correspond to what color
    materialPlottingParameters = Dict{String,PlottingParameters}([
        (
            "Air",
            PlottingParameters(;
                color = [0.95, 0.95, 0.95],
                alpha = 0,
                shade = false,
                lineWidth = 0,
                lineStyle = "None",
            ),
        ),
        (
            "Ag",
            PlottingParameters(;
                color = [0.4, 0.4, 0.4],
                alpha = 0.25,
                shade = true,
                lineStyle = "-",
                lineColor = "r",
            ),
        ),
        (
            "Al2O3",
            PlottingParameters(;
                color = [0.2, 0.2, 0.7],
                alpha = 0.25,
                shade = true,
                lineStyle = "-",
                lineColor = [0.2, 0.7, 0.2],
            ),
        ),
    ])


    # Plot HarmonicsSet:
    # plotHarmonicsSet(simulationDefinition)

    # Plot G-vectorSet:
    # plotGvectorSet(simulationDefinition; scale=μm)

    # Plot x,y components of k-vectors
    # plotkXYVectors(simulationDefinition; scale=μm)

    # Plot lattice:
    # plotLatticeUnit(lattice; scale=μm)

    # Plot reciprocal lattice
    # plotReciprocalLatticeUnit(lattice; scale=μm)

    # Plot the coordinates that are sampled for the convolution matrix
    # plotLayerPositionGrid(layer2, simulationDefinition; scale=μm)

    # Plot the material distribution by color.
    # plotLayerMaterialsDistribution(layer2, simulationDefinition, materialPlottingParameters; scale=μm)
    # plotLayerMaterialsDistribution(layer3, simulationDefinition, materialPlottingParameters; scale=μm)
    # plotLayerMaterialsDistribution(layer4, simulationDefinition, materialPlottingParameters; scale=μm)
    # plotLayerMaterialsDistribution(layer5, simulationDefinition, materialPlottingParameters; scale=μm)
    # plotLayerMaterialsDistribution(layer6, simulationDefinition, materialPlottingParameters; scale=μm)
    # plotLayerMaterialsDistribution(layer7, simulationDefinition, materialPlottingParameters; scale=μm)
    # plotLayerMaterialsDistribution(layer8, simulationDefinition, materialPlottingParameters; scale=μm)


    # Plot x,y components of the zero-order k-vector.
    # plot2DZeroOrderKVector(simulationDefinition; scale=μm)

    # Plot x,y components of the injected k-vectors.
    # # # plot2DinjectedKVectors(simulationDefinition; scale=μm) # Not doing these anymore.  Should just use the E field set separately.


    # Plot 3d k-vectors and polarization vectors
    plot3DinjectedKandPVectors(simulationDefinition; scale=nm, Escale = 500)

    # Plot a 3D structure.
    # plotPatch3D(layerStack, simulationDefinition, materialPlottingParameters; scale=μm)
    # add3DinjectedKandPVectorsToPlot(simulationDefinition; scale=μm, Escale=1)
    # plotPatch3D(layerStack, simulationDefinition, materialPlottingParameters; scale=nm)
    # add3DinjectedKandPVectorsToPlot(simulationDefinition; scale=nm, Escale=500)


    # Plot 2D k-vectors and polarization vectors
    # plot2DinjectedKandPVectors(simulationDefinition; scale=nm, Escale = Escale*500)
    # plot2DinjectedKandPVectors(simulationDefinition; scale=μm, Escale = 0.5)





    # # Plot input and output modes
    middle = UniformLayerDefinition(900 * nm, "Air")
    substrateLayer = SemiInfiniteLayerDefinition("Air")
    superstrateLayer = SemiInfiniteLayerDefinition("Air")
    layerStack = [substrateLayer, middle, superstrateLayer]
    simulationDefinition = SimulationDefinition( lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition )


    # bottomOrders = [_2VectorInt(-2,-2), _2VectorInt(0,0)]
    # topOrders = [_2VectorInt(-2,-2), _2VectorInt(0,0)]
    bottomOrders = [ [0,0], [1,0], [0,1], [-2,-2]]
    topOrders = [ [0,0], [1,0], [0,1], [-2,-2]]
    derivedParameters = DerivedParameters(simulationDefinition)
    allModeData = runSimulation(AllModesAnalysisDefinition(), simulationDefinition)
    # plotPatch3D(layerStack, simulationDefinition, materialPlottingParameters; scale=μm)
    fig, ax = plot3Dlattice(simulationDefinition.lattice, simulationDefinition.layerStack; scale=nm, title="Listed orders")
    add3DlistedKandPVectorsToPlot( allModeData.inputFields, allModeData.outputFields, bottomOrders, topOrders, simulationDefinition, derivedParameters; scale=nm, Escale = 500 )
    xLimits, yLimits = getLatticePlotLimits(simulationDefinition.lattice)
    zLimits = getLayerStackPlotLimits(simulationDefinition.layerStack)
    setCubicAxes(ax, xLimits, yLimits, zLimits; scale=nm)

    # Plot structure with 1D
    # UVstart = [0, 0]
    # UVstop = [1, 1]
    # numDivisions = 20
    # XYstart = convertUVtoXY(lattice, UVstart)
    # XYstop = convertUVtoXY(lattice, UVstop)
    # plotCrossSection(simulationDefinition, XYstart, XYstop, numDivisions, materialPlottingParameters; scale=μm)




    # numDivisionsXY = 20
    # numDivisionsZ = 30
    # positionLineXY = PositionGridXY( CENTERALIGNMENT, XYstart, XYstop, numDivisionsXY)
    # plotCrossSectionNbyArray(simulationDefinition, positionLineXY, numDivisionsZ; scale = μm)
    # plotCrossSectionϵbyArray(simulationDefinition, positionLineXY, numDivisionsZ; scale = μm)


end

end
