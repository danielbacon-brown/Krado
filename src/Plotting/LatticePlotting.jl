
# TODO: Add support for 1D lattice

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
function addLatticeToPlot(lattice::Lattice; scale=1)
    Lx, Ly = calcLatticeBoundaryLine(lattice)
    Lx, Ly = Lx/scale, Ly/scale
    PyPlot.plot(Lx, Ly; color=LATTICECOLOR)
end

# Add to reciprocal lattice plot a set of lines that follow lattice boundaries
function addReciprocalLatticeToPlot(lattice::Lattice; scale=1)
    Gx, Gy = calcReciprocalLatticeBoundaryLine(lattice)
    Gx, Gy = Gx*scale, Gy*scale
    PyPlot.plot(Gx, Gy; color=LATTICECOLOR)
end


function setPlotLimitsAroundLattice(lattice::Lattice, ax; scale=1)
    
    scaleLabel = LENGTHLABEL[scale]
    
    xLimits, yLimits = getLatticePlotLimits(lattice; scale=scale)
    PyPlot.xlim(xLimits[1], xLimits[2])    
    PyPlot.ylim(yLimits[1], yLimits[2])
    
    PyPlot.xlabel("x ($scaleLabel)")
    PyPlot.ylabel("y ($scaleLabel)")    
    
    ax.axis("equal")
    
end

function setPlotLimitsAroundReciprocalLattice(lattice::Lattice, ax; scale=1)
    
    scaleLabel = LENGTHLABEL[scale]
    
    xLimits, yLimits = getReciprocalLatticePlotLimits(lattice; scale=scale)
    PyPlot.xlim(xLimits[1], xLimits[2])    
    PyPlot.ylim(yLimits[1], yLimits[2])
    
    PyPlot.xlabel("x ($(scaleLabel)⁻¹)")
    PyPlot.ylabel("y ($(scaleLabel)⁻¹)")
    
    ax.axis("equal")
    
end


# Plot the real-space periodicity vectors.
function plotLatticeUnit(lattice::Lattice; scale=1)
    
    scaleLabel = LENGTHLABEL[scale]
    
    fig = PyPlot.figure("Lattice. L₁, L₂", figsize=(5,5))
    ax = PyPlot.axes()
    
    addLatticeToPlot(lattice; scale=scale)
        
    setPlotLimitsAroundLattice(lattice, ax; scale=scale)    
end

# Plot the real-space periodicity vectors.
function plotReciprocalLatticeUnit(lattice::Lattice; scale=1)
    
    scaleLabel = LENGTHLABEL[scale]
    
    fig = PyPlot.figure("Reciprocal lattice. G₁, G₂", figsize=(5,5))
    ax = PyPlot.axes()
    
    addReciprocalLatticeToPlot(lattice; scale=scale)    
    setPlotLimitsAroundReciprocalLattice(lattice, ax; scale=scale)

end

# function plot3Dlattice(lattice::Lattice, layerStack::Vector{<:LayerDefinition})
function plot3Dlattice(lattice::Lattice, layerStack::Vector{<:LayerDefinition}; scale=μm)
    
    # Plot 3D lattice:
    xLimits, yLimits = getLatticePlotLimits(lattice)
    zLimits = getLayerStackPlotLimits(layerStack)
    totalThickness = calcTotalThickness(layerStack)

    latticeX, latticeY = calcLatticeBoundaryLine(lattice)
    latticeZbottom = zeros(size(latticeX))
    latticeZtop = ones(size(latticeX))*totalThickness
    latticeX = latticeX/scale
    latticeY = latticeY/scale
    latticeZbottom = latticeZbottom/scale
    latticeZtop = latticeZtop/scale

    # Plot bottom and top:
    plot(latticeX, latticeY, latticeZbottom, color=LATTICECOLOR)
    plot(latticeX, latticeY, latticeZtop, color=LATTICECOLOR)
    # Plot vertical lines:
    for i = 1:4
        plot([latticeX[i],latticeX[i]], [latticeY[i],latticeY[i]], zLimits/scale, color=LATTICECOLOR)
    end
end
