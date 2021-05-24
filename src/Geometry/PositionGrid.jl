

# start represents analytical x,y coordinates of grid start.
# StopU represents analytical x,y coordinates of end of U coordinate
# StopV represents analytical x,y coordinates of end of V coordinate
# StopUV represents analytical x,y coordinates of end of UV coordinate
mutable struct PositionGridXY{N}
    start::_2VectorFloat
    stopU::_2VectorFloat
    stopV::_2VectorFloat
    stopUV::_2VectorFloat
    positions::Array{_2VectorFloat,N}

    function PositionGridXY(start::_2VectorFloat, stopU::_2VectorFloat, stopV::_2VectorFloat, stopUV::_2VectorFloat, positions::Array{_2VectorFloat,N}) where N
        return new{N}(start, stopU, stopV, stopUV, positions)
    end
end

mutable struct PositionGridZ
    start::Float64
    stop::Float64
    positions::Vector{Float64}

    function PositionGridZ(start::Float64, stop::Float64, positions::Vector{Float64})
        return new(start, stop, positions)
    end
end


# Returns a sequence of coordinates at evenly spaced midpoints between xyStart to xyStop.
# Linear
# function PositionGridXYbyMidpoint( start, stop, numDivisions::Integer)
function PositionGridXY( gridAlignment::CenterAlignment, start, stop, numDivisions::Integer)

    start = convert(_2VectorFloat, start)
    stop = convert(_2VectorFloat, stop)

    totalLength = norm(start.-stop)
    pixelLength = totalLength / numDivisions
    # fractionalLength = range0to1exclusiveMidpoint(numDivisions)
    fractionalLength = UnitLinspace(gridAlignment, numDivisions)
    xyPositions = [ start + (stop-start).*fractionalLength[i] for i in 1:numDivisions]
    return PositionGridXY(start, stop, start, stop, xyPositions)
end

# Returns a sequence of coordinates at evenly spaced midpoints between zStart to zStop.
function PositionGridZbyMidpoint( zStart::Real, zStop::Real, numDivisions::Integer)

    totalLength = abs(zStart-zStop)
    pixelLength = totalLength / numDivisions
    # fractionalLength = range0to1exclusiveMidpoint(numDivisions)
    fractionalLength = UnitLinspace(CENTERALIGNMENT, numDivisions)
    zPositions = [ zStart + (zStop-zStart).*fractionalLength[i] for i in 1:numDivisions]
    return PositionGridZ( zStart, zStop, zPositions)
end

# Dispatch over the alignment of the lattice.
function PositionGridXY(lattice::Lattice, numDivisions)
    return PositionGridXY(lattice.gridAlignment, lattice, numDivisions)
end

# 2D lattice
# function PositionGridXYleftAligned( lattice::Lattice, numDivisions::AbstractArray{<:Any})
function PositionGridXY( gridAlignment::LeftAlignment, lattice::Lattice, numDivisions::AbstractArray{<:Any})

    # Fractional distances along U, V vectors
    # Ufractions = range0to1exclusive(numDivisions[U])
    # Vfractions = range0to1exclusive(numDivisions[V])
    Ufractions = unitLinspace(gridAlignment, numDivisions[U])
    Vfractions = unitLinspace(gridAlignment, numDivisions[V])

    gridUVfractions = collect( Base.product(Ufractions, Vfractions) )
    gridUVfractions = reshape(gridUVfractions, Tuple(numDivisions) )

    xyPositions = map( uv -> convertUVtoXY(lattice, _2VectorFloat(uv)), gridUVfractions)

    start = convertUVtoXY(lattice, _2VectorFloat(0,0))
    stopU = convertUVtoXY(lattice, _2VectorFloat(1,0))
    stopV = convertUVtoXY(lattice, _2VectorFloat(0,1))
    stopUV = convertUVtoXY(lattice, _2VectorFloat(1,1))

    return PositionGridXY(start, stopU, stopV, stopUV, xyPositions)
end

function PositionGridXY( gridAlignment::CenterAlignment, lattice::Lattice, numDivisions::AbstractArray{<:Any})

    # Fractional distances along U, V vectors
    # Ufractions = range0to1exclusiveMidpoint(numDivisions[U])
    # Vfractions = range0to1exclusiveMidpoint(numDivisions[V])
    Ufractions = unitLinspace(gridAlignment, numDivisions[U])
    Vfractions = unitLinspace(gridAlignment, numDivisions[V])

    gridUVfractions = collect( Base.product(Ufractions, Vfractions) )
    gridUVfractions = reshape(gridUVfractions, Tuple(numDivisions) )

    xyPositions = map( uv -> convertUVtoXY(lattice, _2VectorFloat(uv)), gridUVfractions)

    start = convertUVtoXY(lattice, _2VectorFloat(0,0))
    stopU = convertUVtoXY(lattice, _2VectorFloat(1,0))
    stopV = convertUVtoXY(lattice, _2VectorFloat(0,1))
    stopUV = convertUVtoXY(lattice, _2VectorFloat(1,1))

    return PositionGridXY(start, stopU, stopV, stopUV, xyPositions)
end

# For 1D lattice
# function PositionGridXYleftAligned( lattice::Lattice, numDivisions::Integer)
function PositionGridXY( gridAlignment::LeftAlignment, lattice::Lattice, numDivisions::Integer)
    @assert(is1D(lattice))
    return PositionGridXY( gridAlignment, lattice, _2VectorInt(numDivisions,1))
end





function Base.length(positionGrid::PositionGridXY)
    return length(positionGrid.positions)
end
function Base.length(positionGrid::PositionGridZ)
    return length(positionGrid.positions)
end

# Only makes sense to get the distance of a linear position grid
function totalDistance(positionGridXY::PositionGridXY{1})
    return norm(positionGridXY.start .- positionGridXY.stopU)
end

function relativeDistances(positionGrid::PositionGridXY)
    return [norm(positionGrid.positions[i] - positionGrid.start) for i in 1:length(positionGrid)]
end
