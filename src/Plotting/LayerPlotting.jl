# Plotting related to layer geometry


function plotLayerPositionGrid(layerDef::PatternedLayerDefinition, simulationDefinition::SimulationDefinition; scale=μm)

    scaleLabel = LENGTHLABEL[scale]

    fig = PyPlot.figure("Layer coordinates", figsize=(5,5))
    ax = PyPlot.axes()
    # PyPlot.clf()

    # Plot lattice unit cell
    addLatticeToPlot(simulationDefinition.lattice; scale=scale)

    # Plot position grid
    # posGrid = calcUniformGridPositions(simulationDefinition.lattice, layerDef)
    posGrid = PositionGridXY(simulationDefinition.lattice, layerDef.numDivisions)
    # posGridValues = posGrid.positions / scale
    xCoords, yCoords = linearizePositionGrid(posGrid)
    # xCoords, yCoords = linearizePositionGrid(posGridValues)

    PyPlot.scatter(xCoords[:]/scale, yCoords[:]/scale, s=1, color="black", marker="," )

    setPlotLimitsAroundLattice(simulationDefinition.lattice, ax; scale=scale)

end


function plotLayerMaterialsDistribution(layerDef::PatternedLayerDefinition, simulationDefinition::SimulationDefinition, materialParams::Dict{String,PlottingParameters}; scale=μm)

    scaleLabel = LENGTHLABEL[scale]

    fig = PyPlot.figure("Layer Materials Distribution", figsize=(5,5))
    ax = PyPlot.axes()

    # Plot lattice unit cell
    # Lx, Ly = calcLatticeBoundaryLine(simulationDefinition.lattice)
    addLatticeToPlot(simulationDefinition.lattice; scale=scale)

    # Get position grid
    # posGrid = calcUniformGridPositions(simulationDefinition.lattice, layerDef)
    posGrid = PositionGridXY(simulationDefinition.lattice, layerDef.numDivisions)
    xCoords, yCoords = linearizePositionGrid(posGrid)
    xCoords = xCoords/scale
    yCoords = yCoords/scale

    materialNameGrid = getMaterialAtPosition( layerDef.layerPattern, posGrid )
    linearMaterialNames = linearizeMaterialNameGrid(materialNameGrid)

    colorGrid = map(materialName -> materialParams[materialName].color, linearMaterialNames)

    PyPlot.scatter(xCoords, yCoords, s=1, c=colorGrid, marker="," )


    setPlotLimitsAroundLattice(simulationDefinition.lattice, ax; scale=scale)

    addMaterialLegend(materialParams::Dict{String,PlottingParameters})

end

function addMaterialLegend(materialParams)
    legendPatches = []
    for (matName, matParam) in materialParams
        push!(legendPatches, PATCHES.Patch(color=matParam.color, label=matName))
    end
    PyPlot.legend(handles=legendPatches, bbox_to_anchor=(1, 1), loc="best")
end
