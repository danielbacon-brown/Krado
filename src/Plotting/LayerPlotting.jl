export plotLayerPositionGrid, plotLayerMaterialsDistribution
# Plotting related to layer geometry


function plotLayerPositionGrid(layerDef::PatternedLayerDefinition, lattice::Lattice; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    fig, ax = create2Dfigure(title="Positions within layer")

    # Plot lattice unit cell
    addLatticeToPlot(ax, lattice; scale=scale)

    # Plot position grid
    posGrid = PositionGridXY(lattice, layerDef.numDivisions)
    xCoords, yCoords = linearizePositionGrid(posGrid)

    ax.scatter(xCoords[:]/scale, yCoords[:]/scale, s=1, color="black", marker="," )

    setPlotLimitsAroundLattice(ax, lattice; scale=scale)

    return fig, ax
end


# function plotLayerMaterialsDistribution(layerDef::PatternedLayerDefinition, lattice::Lattice, materialParams::Dict{String,PlottingParameters}; scale=1)
function plotLayerMaterialsDistribution(layerDef::PatternedLayerDefinition, lattice::Lattice, materialParams::PlottingParameterCollection; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    fig, ax = create2Dfigure(title="Layer Materials Distribution")

    # Plot lattice unit cell
    addLatticeToPlot(ax, lattice; scale=scale)

    # Get position grid
    posGrid = PositionGridXY(lattice, layerDef.numDivisions)
    xCoords, yCoords = linearizePositionGrid(posGrid)
    xCoords = xCoords/scale
    yCoords = yCoords/scale

    materialNameGrid = getMaterialAtPosition( layerDef.layerPattern, posGrid )
    linearMaterialNames = linearizeMaterialNameGrid(materialNameGrid)

    # colorGrid = map(materialName -> materialParams[materialName].color, linearMaterialNames)
    colorGrid = map(materialName -> materialParams[materialName].color, linearMaterialNames)

    ax.scatter(xCoords, yCoords, s=1, c=colorGrid, marker="," )

    setPlotLimitsAroundLattice(ax, lattice; scale=scale)

    # addMaterialLegend(ax, materialParams::Dict{String,PlottingParameters})
    addMaterialLegend(ax, materialParams)

    return fig, ax
end
