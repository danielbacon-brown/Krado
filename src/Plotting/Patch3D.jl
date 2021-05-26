
# function set3DaxisLabels(ax; scale=1)
#     scaleLabel = LENGTHLABEL[scale]
#     ax.xlabel("X ($scaleLabel)")
#     ax.ylabel("Y ($scaleLabel)")
#     ax.zlabel("Z ($scaleLabel)")
# end

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
    # set3DplotLimits(ax, simulationDefinition::SimulationDefinition; scale=scale)
    # setCubicAxes(ax, xLimits, yLimits, zLimits; scale=scale)



    # PyPlot.title("Surface Plot")

    # addMaterialLegend(ax, materialParams::Dict{String,PlottingParameters})

end

function plot3Dstructure(lattice::Lattice, layerStack::Vector{<:LayerDefinition}, materialParams::Dict{String, PlottingParameters}; scale=1)

    # scaleLabel = LENGTHLABEL[scale]

    # fig = figure("3D Materials Distribution", figsize=(5,5))
    # ax = Axes3D(fig)
    #
    # totalThickness = calcTotalThickness(layerStack)

    fig, ax = plot3Dlattice(lattice, layerStack; scale=scale, title="3D Materials Distribution")

    # Plot 3D lattice:
    # xLimits, yLimits = getLatticePlotLimits(simulationDefinition.lattice)
    # zLimits = getLayerStackPlotLimits(layerStack)


    # # Plot substrate:
    # plot3Dsubstrate(ax, lattice, layerStack, materialParams; scale=scale)
    # plot3Dsuperstrate(ax, lattice, layerStack, materialParams; scale=scale)
    # # Plot each layer
    # layerZpositions = calcLayerZpositions(layerStack)
    # for iLayer = 2:(length(layerStack)-1)
    #     plot3Dlayer(ax, lattice, layerStack[iLayer], layerZpositions[iLayer], materialParams; scale=scale)
    # end
    # set3DplotLimits(ax, lattice, layerStack; scale=scale)
    # # set3DplotLimits(ax, simulationDefinition::SimulationDefinition; scale=scale)
    # # setCubicAxes(ax, xLimits, yLimits, zLimits; scale=scale)
    #
    # xlabel("X ($scaleLabel)")
    # ylabel("Y ($scaleLabel)")
    # zlabel("Z ($scaleLabel)")
    #
    # # PyPlot.title("Surface Plot")
    #
    # addMaterialLegend(ax, materialParams::Dict{String,PlottingParameters})

    addToPlot3Dstructure( ax, lattice::Lattice, layerStack::Vector{<:LayerDefinition}, materialParams::Dict{String, PlottingParameters}; scale=scale)

    set3DplotLimits(ax, lattice, layerStack; scale=scale)
    # set3DaxisLabels(ax; scale=scale)

    return fig, ax
end
# plotPatch3D(simulationDefinition::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=1) = plotPatch3D(simulationDefinition.layerStack, simulationDefinition, materialParams; scale=scale)

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
    # @assert numLayers >= 2
    # @assert numLayers >= 3
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

function setCubicAxes(ax, xLimits::Vector{<:Real}, yLimits::Vector{<:Real}, zLimits::Vector{<:Real}; scale=1)
    lengthX = abs(xLimits[2]-xLimits[1])
    lengthY = abs(yLimits[2]-yLimits[1])
    lengthZ = abs(zLimits[2]-zLimits[1])
    centerX = (xLimits[1]+xLimits[2])/2
    centerY = (yLimits[1]+yLimits[2])/2
    centerZ = (zLimits[1]+zLimits[2])/2
    maxLength = max( lengthX, lengthY, lengthZ )
    ax.set_xlim3d( (centerX-maxLength/2)/scale, (centerX+maxLength/2)/scale)
    ax.set_ylim3d( (centerY-maxLength/2)/scale, (centerY+maxLength/2)/scale)
    ax.set_zlim3d( (centerZ-maxLength/2)/scale, (centerZ+maxLength/2)/scale)


    scaleLabel = LENGTHLABEL[scale]
    xlabel("X ($scaleLabel)")
    ylabel("Y ($scaleLabel)")
    zlabel("Z ($scaleLabel)")

end

# function set3DplotLimits(ax, simulationDefinition::SimulationDefinition; scale=1)
function set3DplotLimits(ax, lattice::Lattice, layerStack::Vector{<:LayerDefinition}; scale=1)
    xLimits, yLimits = getLatticePlotLimits(lattice)
    zLimits = getLayerStackPlotLimits(layerStack)
    setCubicAxes(ax, xLimits, yLimits, zLimits; scale=scale)
end


function setCrossSectionAxes(ax, xLimits::Vector{<:Real}, zLimits::Vector{<:Real}; scale=1)
    lengthX = abs(xLimits[2]-xLimits[1])
    lengthZ = abs(zLimits[2]-zLimits[1])
    centerX = (xLimits[1]+xLimits[2])/2
    centerZ = (zLimits[1]+zLimits[2])/2
    maxLength = max( lengthX, lengthZ )
    ax.set_xlim( (centerX-lengthX/2)/scale, (centerX+lengthX/2)/scale)
    ax.set_ylim( (centerZ-lengthZ/2)/scale, (centerZ+lengthZ/2)/scale)

    scaleLabel = LENGTHLABEL[scale]
    xlabel("Distance in X-Y plane ($scaleLabel)")
    ylabel("Z ($scaleLabel)")

end




# 2D cross-section showing material names
# function plotCrossSection(simulationDefinition::SimulationDefinition, XYstart::TU2VectorReal, XYstop::TU2VectorReal, numDivisions::Integer, materialParams::Dict{String, PlottingParameters}; scale=1)
function plotCrossSection(simulationDefinition::SimulationDefinition, positionLine::PositionGridXY, materialParams::Dict{String, PlottingParameters}; scale=1)

    scaleLabel = LENGTHLABEL[scale]


    fig = PyPlot.figure("2D Cross-section Material Distribution", figsize=(5,5))
    ax = PyPlot.axes()

    totalThickness = calcTotalThickness(simulationDefinition.layerStack)

    zLimits = getLayerStackPlotLimits(simulationDefinition.layerStack)
    # linearXYdistance = norm(XYstop - XYstart)
    linearXYdistance = norm(positionLine.stopU - positionLine.start)
    xLimits = [0, linearXYdistance]

    # plotCrossSectionSubstrate(ax, simulationDefinition.layerStack[1], XYstart, XYstop, simulationDefinition, materialParams; scale=scale)
    # plotCrossSectionSuperstrate(ax, last(simulationDefinition.layerStack), XYstart, XYstop, simulationDefinition, materialParams; scale=scale)
    plotCrossSectionSubstrate(ax, simulationDefinition.layerStack[1], xLimits, simulationDefinition, materialParams; scale=scale)
    plotCrossSectionSuperstrate(ax, last(simulationDefinition.layerStack), xLimits, simulationDefinition, materialParams; scale=scale)

    # Plot each layer
    layerZpositions = calcLayerZpositions(simulationDefinition.layerStack)
    for iLayer = 2:(length(simulationDefinition.layerStack)-1)
        # plotCrossSectionLayer(ax, simulationDefinition.layerStack[iLayer], XYstart, XYstop, numDivisions, layerZpositions[iLayer], simulationDefinition, materialParams; scale=scale)
        plotCrossSectionLayer(ax, simulationDefinition.layerStack[iLayer], positionLine, layerZpositions[iLayer], simulationDefinition, materialParams; scale=scale)
    end

    setCrossSectionAxes(ax, xLimits, zLimits; scale=scale)

    # PyPlot.title("Cross-section")

    addMaterialLegend(ax, materialParams::Dict{String,PlottingParameters})

    return fig, ax
end
# plotCrossSection(simulationDefinition::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=1) = plotCrossSection(simulationDefinition.layerStack, simulationDefinition, materialParams; scale=scale)


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

    fig, ax = plotCrossSection(simulationDefinition, XYstart, XYstop, numDivisions, materialPlottingParameters; scale=scale)
    return fig, ax
end



function plotCrossSectionSubstrate(ax, layer::SemiInfiniteLayerDefinition, xLimits::Vector{<:Real}, simulationDefinition::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=1)

    zLimits = getLayerStackPlotLimits(simulationDefinition.layerStack)
    zOuter = zLimits[BOTTOMINDEX]
    zInner = 0.0

    # linearXYdistance = norm(positionLine.stopU - positionLine.Start)
    # linearDistance = linearDistance(positionGridXY)
    linearDist = xLimits[2]-xLimits[1]
    # xLimits = [0, linearXYdistance]

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
    # linearXYdistance = norm(xyStop - xyStart)
    # xLimits = [0, linearXYdistance]

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

    # linearXYdistance = norm(xyStop - xyStart)
    linearDist = linearDistance(positionLine)
    # linearXYinterval = linearDistance / numDivisions
    linearXYinterval = linearSpacing(positionLine)
    xLimits = [0, linearDistance]

    # numDivisions = size(positionLine.positions,1)
    numDivisions = length(positionLine)

    # positionGrid = PositionGridXYbyMidpoint( xyStart, xyStop, numDivisions)
    # positionGrid = PositionGridXY( CENTERALIGNMENT, xyStart, xyStop, numDivisions)

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

    # params = materialParams[layer.backgroundMaterialName]

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

    # posGrid = calcUniformGridPositionsLineMidpoint( xyStart, xyStop, numDivisions)
    posGrid = PositionGridXY(simulationDefinition.lattice, layerDef.numDivisions)


    xyDistances = [norm(posGrid[i]-xyStart) for i in 1:numDivisions]
    # @show xyDistances
    # materialNameGrid = getMaterialAtPosition( layer.layerPattern, posGrid )

    # Only plot when there is a change in material
    startIndex = 1
    # lastMaterialName = materialNameGrid[1]
    lastN = Nvalues[1]
    for i in 2:numDivisions
        currentN = Nvalues[i]
        # currentMaterialName = materialNameGrid[i]
        if currentN != lastN
        # if currentMaterialName != lastMaterialName # If the material is different, plot last block and set new index and material
            xyStart = xyDistances[startIndex] - 0.5*linearXYinterval
            length = linearXYinterval*(i-startIndex)
            # params = materialParams[lastMaterialName]
            # faceColor = colormap( fractionOfRange(currentN,minMaxN) )
            faceColor = colormap( fractionOfRange(lastN,minMaxN) )
            rect = PATCHES.Rectangle([xyStart,zLower]/scale, length/scale, (zUpper-zLower)/scale, facecolor=faceColor, )
            ax.add_patch(rect)
            startIndex = i
            # lastMaterialName = materialNameGrid[i]
            lastN = Nvalues[i]
        end

    end

    # Plot last:
    i = numDivisions+1
    xyStart = xyDistances[startIndex] - 0.5*linearXYinterval
    length = linearXYinterval*(i-startIndex)
    # params = materialParams[lastMaterialName]
    faceColor = colormap( fractionOfRange(lastN,minMaxN) )

    rect = PATCHES.Rectangle([xyStart,zLower]/scale, length/scale, (zUpper-zLower)/scale, facecolor=faceColor)
    ax.add_patch(rect)

end


# 2D cross-section heatmap of results of input function, (real and complex RI)
function plotCrossSectionFuncByArray( func, simulationDefinition::SimulationDefinition, positionGridXY::PositionGridXY, numDivisionsZ::Integer; scale=1, colormap=DEFAULTSEQUENTIALCOLORMAP, title="")
    wavenumber = getWavenumber(simulationDefinition)

    layerStack = simulationDefinition.layerStack

    scaleLabel = LENGTHLABEL[scale]

    # fig = PyPlot.figure(title, figsize=(5,5))
    # ax = PyPlot.axes()
    fig, ax = create2Dfigure(;title=title)

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

    PyPlot.subplot(121)
    PyPlot.imshow( real.(resultByZ), cmap=colormap, extent=(xLimits[1]/scale, xLimits[2]/scale, zLimits[1]/scale, zLimits[2]/scale))
    scaleLabel = LENGTHLABEL[scale]
    PyPlot.xlabel("Lateral distance (x & y) ($scaleLabel)")
    PyPlot.ylabel("z ($scaleLabel)")
    PyPlot.colorbar()
    PyPlot.title("n")

    PyPlot.subplot(122)
    PyPlot.imshow( imag.(resultByZ), cmap=colormap, extent=(xLimits[1]/scale, xLimits[2]/scale, zLimits[1]/scale, zLimits[2]/scale))
    scaleLabel = LENGTHLABEL[scale]
    PyPlot.xlabel("Lateral distance (x & y) ($scaleLabel)")
    PyPlot.ylabel("z ($scaleLabel)")
    PyPlot.colorbar()
    PyPlot.title("k")

    return fig, ax
end

function plotCrossSectionNbyArray(simulationDefinition::SimulationDefinition, positionLineXY::PositionGridXY, numDivisionsZ::Integer; scale = μm)
    println("function n called")
    wavenumber = getWavenumber(simulationDefinition)
    numDivisionsXY = length(positionLineXY)

    function func(layer, positionLineXY)
        ϵμ = getϵμAtPosition( layer, positionLineXY, simulationDefinition.materialCollection, wavenumber)
        nByLayer = [ convert_ϵμ2n.(ϵμ[iPos][EPSILON], ϵμ[iPos][MU]) for iPos = 1:numDivisionsXY]
        return nByLayer
    end
    title = "2D cross-section n,k"
    fig, ax = plotCrossSectionFuncByArray( func, simulationDefinition, positionLineXY, numDivisionsZ; scale=scale, colormap=DEFAULTSEQUENTIALCOLORMAP, title=title)
    return fig, ax
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
    fig, ax = plotCrossSectionFuncByArray( func, simulationDefinition, positionLineXY, numDivisionsZ; scale=scale, colormap=DEFAULTSEQUENTIALCOLORMAP, title=title)
    return fig, ax
end
