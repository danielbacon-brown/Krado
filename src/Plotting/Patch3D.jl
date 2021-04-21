
function plotPatch3D(layerStack::Vector{<:LayerDefinition}, simulationDef::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=μm) 
    
    scaleLabel = LENGTHLABEL[scale]
    
    fig = figure("3D Materials Distribution", figsize=(5,5))
    ax = Axes3D(fig)

    totalThickness = calcTotalThickness(layerStack)
    
    plot3Dlattice(simulationDef.lattice, layerStack; scale=scale)
    
    # Plot 3D lattice:
    xLimits, yLimits = getLatticePlotLimits(simulationDef.lattice)
    zLimits = getLayerStackPlotLimits(layerStack)

    
    # Plot substrate:
    plotPatch3Dsubstrate(layerStack[1], simulationDef, materialParams; scale=scale)
    plotPatch3Dsuperstrate(layerStack[end], simulationDef, materialParams; scale=scale)
    # Plot each layer
    layerZpositions = calcLayerZpositions(layerStack)
    for iLayer = 2:(length(layerStack)-1)
        plotPatch3D(layerStack[iLayer], simulationDef, layerZpositions[iLayer], materialParams; scale=scale)
    end
    
    setCubicAxes(ax, xLimits, yLimits, zLimits; scale=scale)

    xlabel("X ($scaleLabel)")
    ylabel("Y ($scaleLabel)")
    zlabel("Z ($scaleLabel)")

    PyPlot.title("Surface Plot")
    
    addMaterialLegend(materialParams::Dict{String,PlottingParameters}) 
        
end
plotPatch3D(simulationDef::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=μm) = plotPatch3D(simulationDef.layerStack, simulationDef, materialParams; scale=scale) 


#Get z-coordinates of each layer in stack.  Layer 1 (the substrate) 
# function calcLayerZpositions(layerStack::Vector{<:LayerDefinition})
#     layerZ = zeros(size(layerStack))
#     for iLayer = 3:(length(layerStack))
#         layerZ[iLayer] = layerZ[iLayer-1] + layerStack[iLayer-1].thickness
#     end
#     return layerZ
# end

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
    @assert numLayers >= 3
    return sum( [layerStack[i].thickness for i in 2:(numLayers-1)] )
end

function getLayerStackPlotLimits(layerStack::Vector{<:LayerDefinition}; relativePadding=0.5) 
    totalThickness = calcTotalThickness(layerStack)
    padding = totalThickness*relativePadding
    return [-padding, totalThickness+padding]
end

function setCubicAxes(ax, xLimits::Vector{<:Real}, yLimits::Vector{<:Real}, zLimits::Vector{<:Real}; scale=μm)
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

function setCrossSectionAxes(ax, xLimits::Vector{<:Real}, zLimits::Vector{<:Real}; scale=μm)
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
function plotCrossSection(simulationDef::SimulationDefinition, XYstart::TU2VectorReal, XYstop::TU2VectorReal, numDivisions::Integer, materialParams::Dict{String, PlottingParameters}; scale=μm) 
    
    scaleLabel = LENGTHLABEL[scale]
    
    
    fig = PyPlot.figure("2D Cross-section Material Distribution", figsize=(5,5))
    ax = PyPlot.axes()

    totalThickness = calcTotalThickness(simulationDef.layerStack)

    zLimits = getLayerStackPlotLimits(simulationDef.layerStack)
    linearXYdistance = norm(XYstop - XYstart)
    xLimits = [0, linearXYdistance]
    
    plotCrossSectionSubstrate(ax, simulationDef.layerStack[1], XYstart, XYstop, simulationDef, materialParams; scale=scale)
    plotCrossSectionSuperstrate(ax, last(simulationDef.layerStack), XYstart, XYstop, simulationDef, materialParams; scale=scale)

    # Plot each layer
    layerZpositions = calcLayerZpositions(simulationDef.layerStack)
    for iLayer = 2:(length(simulationDef.layerStack)-1)
        plotCrossSectionLayer(ax, simulationDef.layerStack[iLayer], XYstart, XYstop, numDivisions, layerZpositions[iLayer], simulationDef, materialParams; scale=scale)
    end
    
    setCrossSectionAxes(ax, xLimits, zLimits; scale=scale)

    PyPlot.title("Cross-section")
    
    addMaterialLegend(materialParams::Dict{String,PlottingParameters}) 
        
end
plotCrossSection(simulationDef::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=μm) = plotCrossSection(simulationDef.layerStack, simulationDef, materialParams; scale=scale) 


function plotCrossSectionSubstrate(ax, layer::SemiInfiniteLayerDefinition, xyStart::TU2VectorReal, xyStop::TU2VectorReal, simulationDef::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=μm)

    zLimits = getLayerStackPlotLimits(simulationDef.layerStack)
    zOuter = zLimits[BOTTOMINDEX]
    zInner = 0.0
    
    linearXYdistance = norm(xyStop - xyStart)
    xLimits = [0, linearXYdistance]
        
    params = materialParams[layer.backgroundMaterialName]
    
    rect = PATCHES.Rectangle([0,zOuter]/scale, linearXYdistance/scale, (zInner-zOuter)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)

    ax.add_patch(rect)
end

function plotCrossSectionSuperstrate(ax, layer::SemiInfiniteLayerDefinition, xyStart::TU2VectorReal, xyStop::TU2VectorReal, simulationDef::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=μm)

    zLimits = getLayerStackPlotLimits(simulationDef.layerStack)
    zOuter = zLimits[TOPINDEX]
    totalThickness = calcTotalThickness(simulationDef.layerStack)
    zInner = totalThickness
    
    linearXYdistance = norm(xyStop - xyStart)
    xLimits = [0, linearXYdistance]
        
    params = materialParams[layer.backgroundMaterialName]
    
    rect = PATCHES.Rectangle([0,zInner]/scale, linearXYdistance/scale, (zOuter-zInner)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)

    ax.add_patch(rect)
end

# Uniform layer
function plotCrossSectionLayer(ax, layer::UniformLayerDefinition, xyStart::TU2VectorReal, xyStop::TU2VectorReal, numDivisions::Integer, zPosition::Real, simulationDef::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=μm)

    zLower = zPosition
    zUpper = (zPosition+layer.thickness)
        
    linearXYdistance = norm(xyStop - xyStart)
    xLimits = [0, linearXYdistance]
        
    params = materialParams[layer.backgroundMaterialName]
    
    rect = PATCHES.Rectangle([0,zLower]/scale, linearXYdistance/scale, (zUpper-zLower)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)

    ax.add_patch(rect)
end

# function relativeDistance(positionGrid::PositionGridXY) 
#     return norm.(positionGrid .- positionGrid.start)   
# end
# Patterned layer
function plotCrossSectionLayer(ax, layer::PatternedLayerDefinition, xyStart::TU2VectorReal, xyStop::TU2VectorReal, numDivisions::Integer, zPosition::Real, simulationDef::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=μm)

    zLower = zPosition
    zUpper = (zPosition+layer.thickness)
        
    linearXYdistance = norm(xyStop - xyStart)
    linearXYinterval = linearXYdistance / numDivisions
    xLimits = [0, linearXYdistance]
            
    # posGrid = calcUniformGridPositionsLineMidpoint( xyStart, xyStop, numDivisions)
    positionGrid = PositionGridXYbyMidpoint( xyStart, xyStop, numDivisions)

    # xyDistances = [norm(posGrid[i]-xyStart) for i in 1:numDivisions]
    xyDistances = relativeDistances(positionGrid)
    # @show xyDistances
    materialNameGrid = getMaterialAtPosition( layer.layerPattern, positionGrid )
    
    
    
    # Individually plot each pixel
    # for i in 1:numDivisions
    #     xy = xyDistances[i] - 0.5*linearXYinterval
    #     params = materialParams[materialNameGrid[i]]
    #     rect = PATCHES.Rectangle([xy,zLower]/scale, linearXYinterval/scale, (zUpper-zLower)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)
    #     ax.add_patch(rect)
    # end
    
    
    # Only plot when there is a change in material
    startIndex = 1
    lastMaterialName = materialNameGrid[1]
    for i in 2:numDivisions
        currentMaterialName = materialNameGrid[i]
        if currentMaterialName != lastMaterialName # If the material is different, plot last block and set new index and material
            # xyStart = 0 
            blockStart = xyDistances[startIndex] - 0.5*linearXYinterval
            # xyStart = positionGrid.positions[startIndex] - 0.5*linearXYinterval
            length = linearXYinterval*(i-startIndex)
            params = materialParams[lastMaterialName]
            rect = PATCHES.Rectangle([blockStart,zLower]/scale, length/scale, (zUpper-zLower)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)
            ax.add_patch(rect)
            startIndex = i
            lastMaterialName = materialNameGrid[i]
        end
    
    end
    
    # Plot last:
    i = numDivisions+1
    xyStart = xyDistances[startIndex] - 0.5*linearXYinterval
    length = linearXYinterval*(i-startIndex)
    params = materialParams[lastMaterialName]
    rect = PATCHES.Rectangle([xyStart,zLower]/scale, length/scale, (zUpper-zLower)/scale, facecolor=params.color, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle)
     ax.add_patch(rect)
    
end




function plotCrossSectionNsubstrate(ax, Nvalues::Vector{<:Real}, minMaxN::Vector{<:Real}, xyStart::TU2VectorReal, xyStop::TU2VectorReal, simulationDef::SimulationDefinition; scale=μm, colormap=DEFAULTSEQUENTIALCOLORMAP) 

    # Get Z dimensions
    zLimits = getLayerStackPlotLimits(simulationDef.layerStack)
    zOuter = zLimits[BOTTOMINDEX]
    zInner = 0.0

    # Get X,Y dimensions of rectangle
    linearXYdistance = norm(xyStop - xyStart)
    xLimits = [0, linearXYdistance]
        
    faceColor = colormap( fractionOfRange(first(Nvalues),minMaxN) )
    
    rect = PATCHES.Rectangle([0,zOuter]/scale, linearXYdistance/scale, (zInner-zOuter)/scale, facecolor=faceColor)

    ax.add_patch(rect)
end

function plotCrossSectionNsuperstrate(ax, Nvalues::Vector{<:Real}, minMaxN::Vector{<:Real}, xyStart::TU2VectorReal, xyStop::TU2VectorReal, simulationDef::SimulationDefinition; scale=μm, colormap=DEFAULTSEQUENTIALCOLORMAP)

    # Get Z dimensions
    zLimits = getLayerStackPlotLimits(simulationDef.layerStack)
    zOuter = zLimits[TOPINDEX]
    totalThickness = calcTotalThickness(simulationDef.layerStack)
    zInner = totalThickness

    # Get X,Y dimensions of rectangle
    linearXYdistance = norm(xyStop - xyStart)
    xLimits = [0, linearXYdistance]
        
    faceColor = colormap( fractionOfRange(first(Nvalues),minMaxN) )
    
    rect = PATCHES.Rectangle([0,zOuter]/scale, linearXYdistance/scale, (zInner-zOuter)/scale, facecolor=faceColor)

    ax.add_patch(rect)
end

# Uniform layer
function plotCrossSectionNlayer(ax, layer::UniformLayerDefinition,  Nvalues::Vector{<:Real}, minMaxN::Vector{<:Real}, xyStart::TU2VectorReal,  xyStop::TU2VectorReal, numDivisions::Integer, zPosition::Real, simulationDef::SimulationDefinition; scale=μm, colormap=DEFAULTSEQUENTIALCOLORMAP)

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
function plotCrossSectionNlayer(ax, layer::PatternedLayerDefinition,  Nvalues::Vector{<:Real}, minMaxN::Vector{<:Real}, xyStart::TU2VectorReal, xyStop::TU2VectorReal, numDivisions::Integer, zPosition::Real, simulationDef::SimulationDefinition; scale=μm, colormap=DEFAULTSEQUENTIALCOLORMAP)

    zLower = zPosition
    zUpper = (zPosition+layer.thickness)
        
    linearXYdistance = norm(xyStop - xyStart)
    linearXYinterval = linearXYdistance / numDivisions
    xLimits = [0, linearXYdistance]
            
    posGrid = calcUniformGridPositionsLineMidpoint( xyStart, xyStop, numDivisions)

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




# 2D cross-section heatmap of n, (real RI)
# function plotCrossSectionN(simulationDef::SimulationDefinition, XYstart::TU2VectorReal, XYstop::TU2VectorReal, numDivisions::Integer, wavenumber::Wavenumber; scale=μm, colormap=DEFAULTSEQUENTIALCOLORMAP) 
# 
#     scaleLabel = LENGTHLABEL[scale]
# 
#     fig = PyPlot.figure("2D Cross-section n", figsize=(5,5))
#     ax = PyPlot.axes()
# 
#     totalThickness = calcTotalThickness(simulationDef.layerStack)
# 
#     zLimits = getLayerStackPlotLimits(simulationDef.layerStack)
#     linearXYdistance = norm(XYstop - XYstart)
#     xLimits = [0, linearXYdistance]
# 
#     positionGrid = calcUniformGridPositionsLineMidpoint( XYstart, XYstop, numDivisions)
# 
#     #storage of RI-vector for each layer
#     layerNs = Vector{Vector{Float64}}(undef, length(simulationDef.layerStack))
# 
#     # Fill the layerNs array
#     for iLayer = 1:length(simulationDef.layerStack)
#         ϵμ = getϵμAtPosition( simulationDef.layerStack[iLayer], positionGrid, simulationDef.materialCollection, wavenumber)
#         n = [ convert_ϵμ2n.(ϵμ[iPos][EPSILON], ϵμ[iPos][MU]) for iPos = 1:numDivisions]
#         layerNs[iLayer] = real.(n)
#     end
#     # get max and min values for array
#     minMaxN = [minimum(minimum(layerNs)), maximum(maximum(layerNs)) ]
# 
#     layerZpositions = calcLayerZpositions(simulationDef.layerStack)
# 
# 
#     # Plot the array
#     plotCrossSectionNsubstrate(ax, first(layerNs), minMaxN, XYstart, XYstop, simulationDef; scale=scale, colormap=colormap)
#     plotCrossSectionNsuperstrate(ax, last(layerNs), minMaxN, XYstart, XYstop, simulationDef; scale=scale, colormap=colormap)
# 
#     for iLayer = 2:(length(simulationDef.layerStack)-1)
#         plotCrossSectionNlayer(ax, simulationDef.layerStack[iLayer],  layerNs[iLayer], minMaxN, XYstart, XYstop, numDivisions, layerZpositions[iLayer], simulationDef; scale=scale, colormap=colormap)
#     end
# 
#     setCrossSectionAxes(ax, xLimits, zLimits; scale=scale)
# 
#     PyPlot.title("Cross-section")
# 
# end
# plotCrossSectionN(simulationDef::SimulationDefinition, materialParams::Dict{String, PlottingParameters}; scale=μm) = plotCrossSectionN(simulationDef.layerStack, simulationDef, materialParams; scale=scale) 



# 2D cross-section heatmap of n, (real and complex RI)
# function plotCrossSectionNbyArray(simulationDef::SimulationDefinition, XYstart::TU2VectorReal, XYstop::TU2VectorReal, numDivisionsXY::Integer, numDivisionsZ::Integer, wavenumber::Wavenumber; scale=μm, colormap=DEFAULTSEQUENTIALCOLORMAP) 
# 
#     layerStack = simulationDef.layerStack
# 
#     scaleLabel = LENGTHLABEL[scale]
# 
#     fig = PyPlot.figure("2D Cross-section n", figsize=(5,5))
#     ax = PyPlot.axes()
# 
#     totalThickness = calcTotalThickness(layerStack)
# 
#     zLimits = getLayerStackPlotLimits(layerStack)
#     linearXYdistance = norm(XYstop - XYstart)
#     xLimits = [0, linearXYdistance]
# 
#     positionLineXY = calcUniformGridPositionsLineMidpoint( XYstart, XYstop, numDivisionsXY)
#     positionLineZ = reverse(calcUniformGridPositionsLineMidpoint( first(zLimits), last(zLimits), numDivisionsZ) )  # from top to bottom
# 
#     layerZpositions = calcLayerZpositions(layerStack)
# 
#     # fill array by layer
#     nByLayer = zeros(ComplexF64, (length(layerStack), numDivisionsXY))
#     for iLayer = 1:length(layerStack)
#         ϵμ = getϵμAtPosition( layerStack[iLayer], positionLineXY, simulationDef.materialCollection, wavenumber)
#         nByLayer[iLayer,:] = [ convert_ϵμ2n.(ϵμ[iPos][EPSILON], ϵμ[iPos][MU]) for iPos = 1:numDivisionsXY]
#     end
# 
#     # Fill n array
#     nByZ = zeros(ComplexF64, (numDivisionsZ, numDivisionsXY))
#     for iZ = 1:numDivisionsZ
#         zPos = positionLineZ[iZ]
#         iLayer = getLayerIndexFromZ(layerZpositions, zPos)
#         nByZ[iZ, :] = nByLayer[iLayer,:]
#     end
# 
#     PyPlot.subplot(121)
#     PyPlot.imshow( real.(nByZ), cmap=colormap, extent=(xLimits[1]/scale, xLimits[2]/scale, zLimits[1]/scale, zLimits[2]/scale))
#     scaleLabel = LENGTHLABEL[scale]
#     PyPlot.xlabel("Lateral distance (x & y) ($scaleLabel)")
#     PyPlot.ylabel("z ($scaleLabel)")   
#     PyPlot.colorbar()
#     PyPlot.title("n")
# 
#     PyPlot.subplot(122)
#     PyPlot.imshow( imag.(nByZ), cmap=colormap, extent=(xLimits[1]/scale, xLimits[2]/scale, zLimits[1]/scale, zLimits[2]/scale))
#     scaleLabel = LENGTHLABEL[scale]
#     PyPlot.xlabel("Lateral distance (x & y) ($scaleLabel)")
#     PyPlot.ylabel("z ($scaleLabel)")   
#     PyPlot.colorbar()
#     PyPlot.title("k")
# 
# end


# 2D cross-section heatmap of results of input function, (real and complex RI)
function plotCrossSectionFuncByArray( func, simulationDefinition::SimulationDefinition, positionGridXY::PositionGridXY, numDivisionsZ::Integer; scale=μm, colormap=DEFAULTSEQUENTIALCOLORMAP, title="") 
# function plotCrossSectionFuncByArray( func, simulationDef::SimulationDefinition, XYstart::TU2VectorReal, XYstop::TU2VectorReal, numDivisionsXY::Integer, numDivisionsZ::Integer, wavenumber::Wavenumber; scale=μm, colormap=DEFAULTSEQUENTIALCOLORMAP, title="") 
    wavenumber = getWavenumber(simulationDefinition)
    
    layerStack = simulationDefinition.layerStack
    
    scaleLabel = LENGTHLABEL[scale]
    
    fig = PyPlot.figure(title, figsize=(5,5))
    ax = PyPlot.axes()

    totalThickness = calcTotalThickness(layerStack)

    zLimits = getLayerStackPlotLimits(layerStack)
    # linearXYdistance = norm(XYstop - XYstart)
    # linearXYdistance = norm( last(positionLineXY) - first(positionLineXY))
    linearXYdistance = totalDistance(positionGridXY)
    xLimits = [0, linearXYdistance]
    
    # positionLineXY = calcUniformGridPositionsLineMidpoint( XYstart, XYstop, numDivisionsXY)
    positionGridZ = PositionGridZbyMidpoint( last(zLimits), first(zLimits), numDivisionsZ)  # from top to bottom
    # positionLineZ = reverse(calcUniformGridPositionsLineMidpoint( first(zLimits), last(zLimits), numDivisionsZ) )  # from top to bottom
    
    layerZpositions = calcLayerZpositions(layerStack)
    
    # fill result array by layer
    resultByLayer = zeros(ComplexF64, (length(layerStack), length(positionGridXY)))
    for iLayer = 1:length(layerStack)
        
        resultByLayer[iLayer,:] = func(layerStack[iLayer], positionGridXY)
        # ϵμ = getϵμAtPosition( layerStack[iLayer], positionLineXY, simulationDef.materialCollection, wavenumber)
        # nByLayer[iLayer,:] = [ convert_ϵμ2n.(ϵμ[iPos][EPSILON], ϵμ[iPos][MU]) for iPos = 1:numDivisionsXY]
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
    
end

function plotCrossSectionNbyArray(simulationDefinition::SimulationDefinition, positionLineXY::PositionGridXY, numDivisionsZ::Integer; scale = μm)
# function plotCrossSectionNbyArray(simulationDefinition::SimulationDefinition, positionLineXY::Vector{<:TU2VectorReal}, numDivisionsZ::Integer; scale = μm)
    wavenumber = getWavenumber(simulationDefinition)
    numDivisionsXY = length(positionLineXY)

    function func(layer, positionLineXY)
        ϵμ = getϵμAtPosition( layer, positionLineXY, simulationDefinition.materialCollection, wavenumber)
        nByLayer = [ convert_ϵμ2n.(ϵμ[iPos][EPSILON], ϵμ[iPos][MU]) for iPos = 1:numDivisionsXY]
        return nByLayer
    end
    title = "2D cross-section n,k"
    plotCrossSectionFuncByArray( func, simulationDefinition, positionLineXY, numDivisionsZ; scale=μm, colormap=DEFAULTSEQUENTIALCOLORMAP, title=title) 
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
    plotCrossSectionFuncByArray( func, simulationDefinition, positionLineXY, numDivisionsZ; scale=μm, colormap=DEFAULTSEQUENTIALCOLORMAP, title=title) 
    
end