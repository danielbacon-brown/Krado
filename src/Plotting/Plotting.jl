
mutable struct PlottingParameters
    # faces:
    color  # may be a vector or string
    alpha::Float64
    shade::Bool
    # edges:
    lineColor # may be a vector or string
    lineWidth::Float64
    lineStyle::String
    # vertex:
    marker::String

    function PlottingParameters(; color="gray", alpha=0.8, shade=true, lineColor="black", lineWidth=0.5, lineStyle="-", marker="," )
        return new(
            color,
            alpha,
            shade,
            lineColor,
            lineWidth,
            lineStyle,
            marker)
    end
end





# OLD
# IS THIS USED?
# Plot grid of reals
function plotGrid(aGrid::Array{T,2}, titleStr::String, xLabel::String, yLabel::String) where T<:Real

    pygui(true)

    # # Create figure and set axes
    fig = PyPlot.figure(titleStr, figsize=(5,5))
    ax = PyPlot.axes()
    PyPlot.pcolormesh(transpose(aGrid))
    # ax.axis("equal")
    PyPlot.xlim(0,size(aGrid,X))
    PyPlot.ylim(0,size(aGrid,Y))
    PyPlot.colorbar()
    PyPlot.xlabel(xLabel)
    PyPlot.ylabel(yLabel)

    return fig
end

# Plot grid of complex values
function plotGrid(aGrid::Array{T,2}, realTitleStr::String, imagTitleStr::String, xLabel::String, yLabel::String) where {T<:Complex}

    figReal = plotGrid(real(aGrid), realTitleStr, xLabel, yLabel )
    figImag = plotGrid(imag(aGrid), imagTitleStr, xLabel, yLabel )
    return figReal
end





# Makes a scatterplot for the position grid
function plotPositionScatter(posGrid::Array{T,2}, titleStr::String, xLabel::String, yLabel::String) where T<:TU2VectorReal
    pygui(true)

    fig = PyPlot.figure(titleStr, figsize=(5,5))
    ax = PyPlot.axes()

    xCoords, yCoords = linearizePositionGrid(posGrid)

    PyPlot.scatter(xCoords[:], yCoords[:])
    PyPlot.title(titleStr)
    PyPlot.xlabel("x")
    PyPlot.ylabel("y")
    return fig
end

# Scatterplot of position grid of the given layer
function plotPositionScatter(layerDef::LayerDefinition, lattice::Lattice)
    posGrid = calcUniformGridPositions(lattice, layerDef)
    return plotPositionScatter(posGrid, "Layer Positions", "X", "Y")
end


function plotGridByScatter(aGrid::Array{T1,2}, posGrid::Array{_2VectorFloat,2}, titleStr::String, xLabel::String, yLabel::String) where {T1<:Real}

    pygui(true)
    fig = PyPlot.figure(titleStr, figsize=(5,5))
    ax = PyPlot.axes()

    # xCoords, yCoords = linearizePositionGrid(posGrid)
    xGrid, yGrid = extractPositionGridComponents(posGrid)

    PyPlot.scatter(xGrid, yGrid, c=aGrid )
    PyPlot.colorbar()
    PyPlot.xlabel(xLabel)
    PyPlot.ylabel(yLabel)

end

plotGridByScatter(aGrid::Array{T1,2}, posGrid::Array{T2,2}, titleStr::String, xLabel::String, yLabel::String) where {T1<:Real, T2<:TU2VectorReal} = plotGridByScatter(aGrid, convert(Array{_2VectorFloat,2},posGrid), titleStr, xLabel, yLabel)
