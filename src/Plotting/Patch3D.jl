

function addToPlot3Dstructure( ax, lattice::Lattice, layerStack::Vector{<:LayerDefinition}, materialParams::Dict{String, PlottingParameters}; scale=1)
    # Plot substrate:
    plot3Dsubstrate(ax, lattice, layerStack, materialParams; scale=scale)
    plot3Dsuperstrate(ax, lattice, layerStack, materialParams; scale=scale)
    # Plot each layer
    layerZpositions = calcLayerZpositions(layerStack)
    for iLayer = 2:(length(layerStack)-1)
        plot3Dlayer(ax, lattice, layerStack[iLayer], layerZpositions[iLayer], materialParams; scale=scale)
    end
    set3DplotLimits(ax, lattice, layerStack; scale=scale)

end

function plot3Dstructure(lattice::Lattice, layerStack::Vector{<:LayerDefinition}, materialParams::Dict{String, PlottingParameters}; scale=1)

    fig, ax = plot3Dlattice(lattice, layerStack; scale=scale, title="3D Materials Distribution")

    addToPlot3Dstructure( ax, lattice::Lattice, layerStack::Vector{<:LayerDefinition}, materialParams::Dict{String, PlottingParameters}; scale=scale)

    set3DplotLimits(ax, lattice, layerStack; scale=scale)

    return fig, ax
end

#Get z-coordinates of the bottom of each layer in stack.  Layer 1 (the substrate) is below 0, layer 2 is at zero.  Last layer is total thickness
function calcLayerZpositions(layerStack::Vector{<:LayerDefinition})
    zLimits = getLayerStackPlotLimits(layerStack)

    layerZ = zeros(length(layerStack) + 1)
    for iLayer = 3:(length(layerStack))
        layerZ[iLayer] = layerZ[iLayer-1] + layerStack[iLayer-1].thickness
    end
    layerZ[1] = zLimits[1]
    layerZ[length(layerStack) + 1] = zLimits[2]

    return layerZ
end

function calcTotalThickness(layerStack::Vector{<:LayerDefinition})
    numLayers = length(layerStack)

    if numLayers >= 3

        return sum( [layerStack[i].thickness for i in 2:(numLayers-1)] )
    else
        return 0
    end
end

function getLayerStackPlotLimits(layerStack::Vector{<:LayerDefinition}; relativePadding=0.5)
    totalThickness = calcTotalThickness(layerStack)
    padding = totalThickness*relativePadding
    return [-padding, totalThickness+padding]
end








# 2D cross-section showing material names
function plotCrossSection(simulationDefinition::SimulationDefinition, positionLine::PositionGridXY, materialParams::Dict{String, PlottingParameters}; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    totalThickness = calcTotalThickness(simulationDefinition.layerStack)

    zLimits = getLayerStackPlotLimits(simulationDefinition.layerStack)

    linearXYdistance = norm(positionLine.stopU - positionLine.start)
    xLimits = [0, linearXYdistance]

    fig, ax = plt.subplots(1)
    fig.canvas.set_window_title("Cross-section Material Distribution")
    ax.set_xlim( (xLimits[1]/scale, xLimits[2]/scale) )
    ax.set_ylim( (zLimits[1]/scale, zLimits[2]/scale) )
    ax.set_aspect("equal")

    plotCrossSectionSubstrate(ax, simulationDefinition.layerStack[1], xLimits, simulationDefinition, materialParams; scale=scale)
    plotCrossSectionSuperstrate(ax, last(simulationDefinition.layerStack), xLimits, simulationDefinition, materialParams; scale=scale)

    # Plot each layer
    layerZpositions = calcLayerZpositions(simulationDefinition.layerStack)
    for iLayer = 2:(length(simulationDefinition.layerStack)-1)
        plotCrossSectionLayer(ax, simulationDefinition.layerStack[iLayer], positionLine, layerZpositions[iLayer], simulationDefinition, materialParams; scale=scale)
    end

    addMaterialLegend(ax, materialParams::Dict{String,PlottingParameters})
    return fig, ax
end


# Sugar for plotting 1D cross-section.
function plotCrossSection(simulationDefinition::SimulationDefinition, numDivisions::Integer, materialParams::Dict{String, PlottingParameters}; scale=1)
    if ~is1D(simulationDefinition.lattice)
        error("Expected 1-dimensional lattice.  Start and stop positions must be defined for plotting cross-section of 2D structure.")
    end

    UVstart = [0, 0]
    UVstop = [1, 0]
    numDivisions = 20
    XYstart = convertUVtoXY(lattice, UVstart)
    XYstop = convertUVtoXY(lattice, UVstop)
    positionLineXY = PositionGridXY( CENTERALIGNMENT, XYstart, XYstop, numDivisions)

    fig, ax = plotCrossSection(simulationDefinition, positionLineXY, materialPlottingParameters; scale=scale)
    return fig, ax
end



function plotCrossSectionSubstrate(ax, layer::SemiInfiniteLayerDefinition, xLimits::Vector{<:Real}, simulationDefinition::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=1)

    zLimits = getLayerStackPlotLimits(simulationDefinition.layerStack)
    zOuter = zLimits[BOTTOMINDEX]
    zInner = 0.0

    linearDist = xLimits[2]-xLimits[1]

    params = materialParams[layer.backgroundMaterialName]

    rect = PATCHES.Rectangle([0,zOuter]/scale, linearDist/scale, (zInner-zOuter)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)

    ax.add_patch(rect)
end

function plotCrossSectionSuperstrate(ax, layer::SemiInfiniteLayerDefinition, xLimits::Vector{<:Real}, simulationDefinition::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=1)

    zLimits = getLayerStackPlotLimits(simulationDefinition.layerStack)
    zOuter = zLimits[TOPINDEX]
    totalThickness = calcTotalThickness(simulationDefinition.layerStack)
    zInner = totalThickness

    linearDist = xLimits[2]-xLimits[1]

    params = materialParams[layer.backgroundMaterialName]

    rect = PATCHES.Rectangle([0,zInner]/scale, linearDist/scale, (zOuter-zInner)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)

    ax.add_patch(rect)
end

# Uniform layer
function plotCrossSectionLayer(ax, layer::UniformLayerDefinition, positionLine::PositionGridXY, zPosition::Real, simulationDefinition::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=1)

    zLower = zPosition
    zUpper = (zPosition+layer.thickness)

    linearDist = linearDistance(positionLine)
    xLimits = [0, linearDist]

    params = materialParams[layer.backgroundMaterialName]

    rect = PATCHES.Rectangle([0,zLower]/scale, linearDist/scale, (zUpper-zLower)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)

    ax.add_patch(rect)
end


# Patterned layer
function plotCrossSectionLayer(ax, layer::PatternedLayerDefinition, positionLine::PositionGridXY, zPosition::Real, simulationDefinition::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=1)

    zLower = zPosition
    zUpper = (zPosition+layer.thickness)

    linearDist = linearDistance(positionLine)
    linearXYinterval = linearSpacing(positionLine)
    xLimits = [0, linearDistance]

    numDivisions = length(positionLine)

    xyDistances = relativeDistances(positionLine)
    materialNameGrid = getMaterialAtPosition( layer.layerPattern, positionLine )


    # Only plot when there is a change in material
    startIndex = 1
    lastMaterialName = materialNameGrid[1]
    for i in 2:numDivisions
        currentMaterialName = materialNameGrid[i]
        if currentMaterialName != lastMaterialName # If the material is different, plot last block and set new index and material
            blockStart = xyDistances[startIndex] - 0.5*linearXYinterval
            smallLength = linearXYinterval*(i-startIndex)
            params = materialParams[lastMaterialName]
            rect = PATCHES.Rectangle([blockStart,zLower]/scale, smallLength/scale, (zUpper-zLower)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)
            ax.add_patch(rect)
            startIndex = i
            lastMaterialName = materialNameGrid[i]
        end

    end

    # Plot last:
    i = numDivisions+1
    xyStart = xyDistances[startIndex] - 0.5*linearXYinterval
    smallLength = linearXYinterval*(i-startIndex)
    params = materialParams[lastMaterialName]
    rect = PATCHES.Rectangle([xyStart,zLower]/scale, smallLength/scale, (zUpper-zLower)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)
     ax.add_patch(rect)

end




function plotCrossSectionNsubstrate(ax, Nvalues::Vector{<:Real}, minMaxN::Vector{<:Real}, xyStart::TU2VectorReal, xyStop::TU2VectorReal, simulationDefinition::SimulationDefinition; scale=1, colormap=DEFAULTSEQUENTIALCOLORMAP)

    # Get Z dimensions
    zLimits = getLayerStackPlotLimits(simulationDefinition.layerStack)
    zOuter = zLimits[BOTTOMINDEX]
    zInner = 0.0

    # Get X,Y dimensions of rectangle
    linearXYdistance = norm(xyStop - xyStart)
    xLimits = [0, linearXYdistance]

    faceColor = colormap( fractionOfRange(first(Nvalues),minMaxN) )

    rect = PATCHES.Rectangle([0,zOuter]/scale, linearXYdistance/scale, (zInner-zOuter)/scale, facecolor=faceColor)

    ax.add_patch(rect)
end

function plotCrossSectionNsuperstrate(ax, Nvalues::Vector{<:Real}, minMaxN::Vector{<:Real}, xyStart::TU2VectorReal, xyStop::TU2VectorReal, simulationDefinition::SimulationDefinition; scale=1, colormap=DEFAULTSEQUENTIALCOLORMAP)

    # Get Z dimensions
    zLimits = getLayerStackPlotLimits(simulationDefinition.layerStack)
    zOuter = zLimits[TOPINDEX]
    totalThickness = calcTotalThickness(simulationDefinition.layerStack)
    zInner = totalThickness

    # Get X,Y dimensions of rectangle
    linearXYdistance = norm(xyStop - xyStart)
    xLimits = [0, linearXYdistance]

    faceColor = colormap( fractionOfRange(first(Nvalues),minMaxN) )

    rect = PATCHES.Rectangle([0,zOuter]/scale, linearXYdistance/scale, (zInner-zOuter)/scale, facecolor=faceColor)

    ax.add_patch(rect)
end

# Uniform layer
function plotCrossSectionNlayer(ax, layer::UniformLayerDefinition,  Nvalues::Vector{<:Real}, minMaxN::Vector{<:Real}, xyStart::TU2VectorReal,  xyStop::TU2VectorReal, numDivisions::Integer, zPosition::Real, simulationDefinition::SimulationDefinition; scale=1, colormap=DEFAULTSEQUENTIALCOLORMAP)

    zLower = zPosition
    zUpper = (zPosition+layer.thickness)

    linearXYdistance = norm(xyStop - xyStart)
    xLimits = [0, linearXYdistance]

    faceColor = colormap( fractionOfRange(first(Nvalues),minMaxN) )

    rect = PATCHES.Rectangle([0,zLower]/scale, linearXYdistance/scale, (zUpper-zLower)/scale, facecolor=faceColor)

    ax.add_patch(rect)
end

# Patterned layer
function plotCrossSectionNlayer(ax, layer::PatternedLayerDefinition,  Nvalues::Vector{<:Real}, minMaxN::Vector{<:Real}, xyStart::TU2VectorReal, xyStop::TU2VectorReal, numDivisions::Integer, zPosition::Real, simulationDefinition::SimulationDefinition; scale=1, colormap=DEFAULTSEQUENTIALCOLORMAP)

    zLower = zPosition
    zUpper = (zPosition+layer.thickness)

    linearXYdistance = norm(xyStop - xyStart)
    linearXYinterval = linearXYdistance / numDivisions
    xLimits = [0, linearXYdistance]

    posGrid = PositionGridXY(simulationDefinition.lattice, layerDef.numDivisions)


    xyDistances = [norm(posGrid[i]-xyStart) for i in 1:numDivisions]

    # Only plot when there is a change in material
    startIndex = 1
    # lastMaterialName = materialNameGrid[1]
    lastN = Nvalues[1]
    for i in 2:numDivisions
        currentN = Nvalues[i]
        if currentN != lastN
            xyStart = xyDistances[startIndex] - 0.5*linearXYinterval
            length = linearXYinterval*(i-startIndex)
            faceColor = colormap( fractionOfRange(lastN,minMaxN) )
            rect = PATCHES.Rectangle([xyStart,zLower]/scale, length/scale, (zUpper-zLower)/scale, facecolor=faceColor, )
            ax.add_patch(rect)
            startIndex = i
            lastN = Nvalues[i]
        end

    end

    # Plot last:
    i = numDivisions+1
    xyStart = xyDistances[startIndex] - 0.5*linearXYinterval
    length = linearXYinterval*(i-startIndex)
    faceColor = colormap( fractionOfRange(lastN,minMaxN) )

    rect = PATCHES.Rectangle([xyStart,zLower]/scale, length/scale, (zUpper-zLower)/scale, facecolor=faceColor)
    ax.add_patch(rect)

end


# 2D cross-section heatmap of results of input function, (real and complex RI)
function plotCrossSectionFuncByArray( func, simulationDefinition::SimulationDefinition, positionGridXY::PositionGridXY, numDivisionsZ::Integer; scale=1, colormap=DEFAULTSEQUENTIALCOLORMAP, title="")
    wavenumber = getWavenumber(simulationDefinition)

    layerStack = simulationDefinition.layerStack

    scaleLabel = LENGTHLABEL[scale]

    totalThickness = calcTotalThickness(layerStack)

    zLimits = getLayerStackPlotLimits(layerStack)

    linearXYdistance = linearDistance(positionGridXY)
    xLimits = [0, linearXYdistance]

    positionGridZ = PositionGridZbyMidpoint( last(zLimits), first(zLimits), numDivisionsZ)  # from top to bottom

    layerZpositions = calcLayerZpositions(layerStack)

    # fill result array by layer
    resultByLayer = zeros(ComplexF64, (length(layerStack), length(positionGridXY)))
    for iLayer = 1:length(layerStack)
        resultByLayer[iLayer,:] = func(layerStack[iLayer], positionGridXY)
    end

    # Fill result array by z
    resultByZ = zeros(ComplexF64, (numDivisionsZ, length(positionGridXY)))
    for iZ = 1:length(positionGridZ)
        zPos = positionGridZ.positions[iZ]
        iLayer = getLayerIndexFromZ(layerZpositions, zPos)
        resultByZ[iZ, :] = resultByLayer[iLayer,:]
    end

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.canvas.set_window_title(title)

    colormap1 = ax1.imshow( real.(resultByZ), cmap=colormap, extent=(xLimits[1]/scale, xLimits[2]/scale, zLimits[1]/scale, zLimits[2]/scale))

    setLabels(ax1, "Lateral distance (x & y) ($scaleLabel)", "z ($scaleLabel)")

    fig.colorbar(colormap1, ax=ax1)
    ax1.set_title("n")

    colormap2 = ax2.imshow( imag.(resultByZ), cmap=colormap, extent=(xLimits[1]/scale, xLimits[2]/scale, zLimits[1]/scale, zLimits[2]/scale))

    setLabels(ax2, "Lateral distance (x & y) ($scaleLabel)", "z ($scaleLabel)")
    fig.colorbar(colormap2, ax=ax2)
    ax2.set_title("k")

    return fig, ax1, ax2
end

function plotCrossSectionNbyArray(simulationDefinition::SimulationDefinition, positionLineXY::PositionGridXY, numDivisionsZ::Integer; scale = μm)
    wavenumber = getWavenumber(simulationDefinition)
    numDivisionsXY = length(positionLineXY)

    function func(layer, positionLineXY)
        ϵμ = getϵμAtPosition( layer, positionLineXY, simulationDefinition.materialCollection, wavenumber)
        nByLayer = [ convert_ϵμ2n.(ϵμ[iPos][EPSILON], ϵμ[iPos][MU]) for iPos = 1:numDivisionsXY]
        return nByLayer
    end
    title = "2D cross-section n,k"
    fig, ax1, ax2 = plotCrossSectionFuncByArray( func, simulationDefinition, positionLineXY, numDivisionsZ; scale=scale, colormap=DEFAULTSEQUENTIALCOLORMAP, title=title)
    return fig, ax1, ax2
end

function plotCrossSectionϵbyArray(simulationDefinition::SimulationDefinition, positionLineXY::PositionGridXY, numDivisionsZ::Integer; scale = μm)

    wavenumber = getWavenumber(simulationDefinition)
    numDivisionsXY = length(positionLineXY)

    function func(layer, positionLineXY)
        ϵμ = getϵμAtPosition( layer, positionLineXY, simulationDefinition.materialCollection, wavenumber)
        ϵByLayer = [ ϵμ[iPos][EPSILON] for iPos = 1:numDivisionsXY]
        return ϵByLayer
    end
    title = "2D cross-section ϵ\',ϵ\'\'"
    fig, ax1, ax2 = plotCrossSectionFuncByArray( func, simulationDefinition, positionLineXY, numDivisionsZ; scale=scale, colormap=DEFAULTSEQUENTIALCOLORMAP, title=title)
    return fig, ax1, ax2
end
