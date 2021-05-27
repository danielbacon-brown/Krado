

function create2Dfigure(;title::String = "")

    fig = PyPlot.figure(title, figsize=DEFAULTPLOTSIZE)
    ax = PyPlot.axes()

    return fig, ax
end

function create3Dfigure(;title::String = "")

    fig = figure(title, figsize=DEFAULTPLOTSIZE)
    ax = Axes3D(fig)

    return fig, ax
end

function set2Dlimits(ax, xLimits, yLimits)
    ax.set_xlim(xLimits[1], xLimits[2])
    ax.set_ylim(yLimits[1], yLimits[2])
    ax.axis("equal")
end

function setLabels(ax, xLabel, yLabel)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
end

function addMaterialLegend(ax, materialParams)
    legendPatches = []
    for (matName, matParam) in materialParams
        push!(legendPatches, PATCHES.Patch(color=matParam.color, label=matName))
    end
    ax.legend(handles=legendPatches, bbox_to_anchor=(1, 1), loc="best")
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

function set3DplotLimits(ax, lattice::Lattice, layerStack::Vector{<:LayerDefinition}; scale=1)
    xLimits, yLimits = getLatticePlotLimits(lattice)
    zLimits = getLayerStackPlotLimits(layerStack)
    setCubicAxes(ax, xLimits, yLimits, zLimits; scale=scale)
end


function setCrossSectionAxes(ax, xLimits::Vector{<:Real}, zLimits::Vector{<:Real}; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    setLabels(ax, "Distance in X-Y plane ($scaleLabel)", "Z ($scaleLabel)")

    set2Dlimits(ax, xLimits, zLimits)

    PyPlot.tight_layout()

end
