# Plotting related to layer geometry


# function plotLayerPositionGrid(layerDef::PatternedLayerDefinition, simulationDefinition::SimulationDefinition; scale=1)
function plotLayerPositionGrid(layerDef::PatternedLayerDefinition, lattice::Lattice; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    # fig = PyPlot.figure("Layer coordinates", figsize=(5,5))
    # ax = PyPlot.axes()
    # PyPlot.clf()
    fig, ax = create2Dfigure(title="Positions within layer")

    # Plot lattice unit cell
    addLatticeToPlot(ax, lattice; scale=scale)

    # Plot position grid
    # posGrid = calcUniformGridPositions(simulationDefinition.lattice, layerDef)
    posGrid = PositionGridXY(lattice, layerDef.numDivisions)
    # posGridValues = posGrid.positions / scale
    xCoords, yCoords = linearizePositionGrid(posGrid)
    # xCoords, yCoords = linearizePositionGrid(posGridValues)

    PyPlot.scatter(xCoords[:]/scale, yCoords[:]/scale, s=1, color="black", marker="," )

    setPlotLimitsAroundLattice(ax, lattice; scale=scale)

    return fig, ax
end


# function plotLayerMaterialsDistribution(layerDef::PatternedLayerDefinition, simulationDefinition::SimulationDefinition, materialParams::Dict{String,PlottingParameters}; scale=1)
function plotLayerMaterialsDistribution(layerDef::PatternedLayerDefinition, lattice::Lattice, materialParams::Dict{String,PlottingParameters}; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    # fig = PyPlot.figure("Layer Materials Distribution", figsize=(5,5))
    # ax = PyPlot.axes()
    fig, ax = create2Dfigure(title="Layer Materials Distribution")

    # Plot lattice unit cell
    # Lx, Ly = calcLatticeBoundaryLine(simulationDefinition.lattice)
    addLatticeToPlot(ax, lattice; scale=scale)

    # Get position grid
    # posGrid = calcUniformGridPositions(simulationDefinition.lattice, layerDef)
    posGrid = PositionGridXY(lattice, layerDef.numDivisions)
    xCoords, yCoords = linearizePositionGrid(posGrid)
    xCoords = xCoords/scale
    yCoords = yCoords/scale

    materialNameGrid = getMaterialAtPosition( layerDef.layerPattern, posGrid )
    linearMaterialNames = linearizeMaterialNameGrid(materialNameGrid)

    colorGrid = map(materialName -> materialParams[materialName].color, linearMaterialNames)

    PyPlot.scatter(xCoords, yCoords, s=1, c=colorGrid, marker="," )


    setPlotLimitsAroundLattice(ax, lattice; scale=scale)

    addMaterialLegend(ax, materialParams::Dict{String,PlottingParameters})

    return fig, ax
end

function addMaterialLegend(ax, materialParams)
    legendPatches = []
    for (matName, matParam) in materialParams
        push!(legendPatches, PATCHES.Patch(color=matParam.color, label=matName))
    end
    # PyPlot.legend(handles=legendPatches, bbox_to_anchor=(1, 1), loc="best")
    ax.legend(handles=legendPatches, bbox_to_anchor=(1, 1), loc="best")
end
