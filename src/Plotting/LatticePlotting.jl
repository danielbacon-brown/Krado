
# TODO: Add support for 1D lattice

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

# Returns the X- and Y- limits for a plot that's slightly larger than the lattice.
function getLatticePlotLimits(lattice::Lattice; relativePadding=0.1, scale=1)
    Lx, Ly = calcLatticeBoundaryLine(lattice)
    Lx, Ly = Lx/scale, Ly/scale
    paddingX = relativePadding*(maximum(Lx) - minimum(Lx))
    paddingY = relativePadding*(maximum(Ly) - minimum(Ly))
    xLimits = [minimum(Lx)-paddingX, maximum(Lx)+paddingX]
    yLimits = [minimum(Ly)-paddingY, maximum(Ly)+paddingY]
    return xLimits, yLimits
end

# Returns the X- and Y- limits for a plot that's slightly larger than the reciprocal lattice.
function getReciprocalLatticePlotLimits(lattice::Lattice; relativePadding=0.1, scale=1)
    Gx, Gy = calcReciprocalLatticeBoundaryLine(lattice)
    Gx, Gy = Gx*scale, Gy*scale
    paddingX = relativePadding*(maximum(Gx) - minimum(Gx))
    paddingY = relativePadding*(maximum(Gy) - minimum(Gy))
    xLimits = [minimum(Gx)-paddingX, maximum(Gx)+paddingX]
    yLimits = [minimum(Gy)-paddingY, maximum(Gy)+paddingY]
    return xLimits, yLimits
end


# Add to plot a set of lines that follow lattice boundaries
function addLatticeToPlot(ax, lattice::Lattice; scale=1)
    Lx, Ly = calcLatticeBoundaryLine(lattice)
    Lx, Ly = Lx/scale, Ly/scale
    ax.plot(Lx, Ly; color=LATTICECOLOR)
end

# Add to reciprocal lattice plot a set of lines that follow lattice boundaries
function addReciprocalLatticeToPlot(ax, lattice::Lattice; scale=1)
    Gx, Gy = calcReciprocalLatticeBoundaryLine(lattice)
    Gx, Gy = Gx*scale, Gy*scale
    ax.plot(Gx, Gy; color=LATTICECOLOR)
end


function setPlotLimitsAroundLattice(ax, lattice::Lattice; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    xLimits, yLimits = getLatticePlotLimits(lattice; scale=scale)
    set2Dlimits(ax, xLimits, yLimits)
    setLabels(ax, "x ($(scaleLabel))", "y ($(scaleLabel))")
    ax.axis("equal")

end

function setPlotLimitsAroundReciprocalLattice(ax, lattice::Lattice; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    xLimits, yLimits = getReciprocalLatticePlotLimits(lattice; scale=scale)
    set2Dlimits(ax, xLimits, yLimits)
    setLabels(ax, "x ($(scaleLabel)⁻¹)", "y ($(scaleLabel)⁻¹)")
    ax.axis("equal")

end


# Plot the real-space periodicity vectors.
function plotLatticeUnit(lattice::Lattice; scale=1)

    scaleLabel = LENGTHLABEL[scale]
    fig, ax = create2Dfigure(title="Lattice. L₁, L₂")
    addLatticeToPlot(ax, lattice; scale=scale)
    setPlotLimitsAroundLattice(ax, lattice; scale=scale)

    return fig, ax
end

# Plot the real-space periodicity vectors.
function plotReciprocalLatticeUnit(lattice::Lattice; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    fig, ax = create2Dfigure(title="Reciprocal lattice. G₁, G₂")
    addReciprocalLatticeToPlot(ax, lattice; scale=scale)
    setPlotLimitsAroundReciprocalLattice(ax, lattice; scale=scale)

    return fig, ax
end

function addToPlot3Dlattice(ax, lattice::Lattice, layerStack::Vector{<:LayerDefinition}; scale=1)

    # Get boundary limits of plot
    xLimits, yLimits = getLatticePlotLimits(lattice)
    zLimits = getLayerStackPlotLimits(layerStack)
    totalThickness = calcTotalThickness(layerStack)

    # Calcuulate coordinates of lattice
    latticeX, latticeY = calcLatticeBoundaryLine(lattice)
    latticeZbottom = zeros(size(latticeX))
    latticeZtop = ones(size(latticeX))*totalThickness

    # Plot horizontal lines
    plot(latticeX./scale, latticeY./scale, latticeZbottom./scale, color=LATTICECOLOR)
    plot(latticeX./scale, latticeY./scale, latticeZtop./scale, color=LATTICECOLOR)

    # Plot vertical lines
    for i = 1:4
        plot([latticeX[i],latticeX[i]]./scale, [latticeY[i],latticeY[i]]./scale, zLimits./scale, color=LATTICECOLOR)
    end
end


# function plot3Dlattice(lattice::Lattice, layerStack::Vector{<:LayerDefinition})
function plot3Dlattice(lattice::Lattice, layerStack::Vector{<:LayerDefinition}; scale=1, title="Lattice")

    fig, ax = create3Dfigure(title=title)

    addToPlot3Dlattice(ax, lattice::Lattice, layerStack::Vector{<:LayerDefinition}; scale=scale)

    set3DplotLimits(ax, lattice, layerStack; scale=scale)

    return fig, ax
end
