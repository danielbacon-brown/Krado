# Definition of lattice/period

export Lattice


abstract type AbstractLattice end

# Vectors defining the real space and reciprocal space periodicity vectors [x,y] defining the lattice.
# L is real space.  G is reciprocal space.
# originOffset is a real-space shift of the unit cell origin, in terms of [fraction of L1, fraction of L2].  E.g. for center of unit cell. originOffset = [-0.5,-0.5]


"""
    Lattice

Describes the periodicity of the lattice and the position of the unit cell.

L₁ and L₂: real-space vectors describing the periodicity of the lattice.  Uses the X,Y coordinate system.  For 1D periodicity, L₁ contains the periodicity vector and L₂ is 0.
    The periodicity of the lattice is given by

G₁ and G₂: reciprocal space lattice vectors.  Uses reciprocal X,Y coordinates.

originOffsetUV: Vector describing how the unit cell should be translated.  UV coordinates.  Defaults to zero, unless originOffsetXY has been set.
originOffsetXY: Vector describing how the unit cell should be translated.  XY coordinates.  Defaults to zero, unless originOffsetUV has been set.
If only one originOffset parameter has been set, that values is used to calculate the other one.  If neither originOffset parameter has been set, both are [0,0].
"""
mutable struct Lattice <: AbstractLattice
    # Real space lattice vectors
    L₁::_2VectorFloat
    L₂::_2VectorFloat

    # Reciprocal space lattice vectors
    G₁::_2VectorFloat
    G₂::_2VectorFloat

    # Location of the origin of the unit cell (in UV coordinates)
    # originOffset::_2VectorFloat
    originOffsetUV::_2VectorFloat
    originOffsetXY::_2VectorFloat

    # Transformation matrices to convert between XY and UV coordinates.  (Real space)
    transformationUVtoXY::Array{Float64,2}
    transformationXYtoUV::Array{Float64,2}

    gridAlignment::GridAlignment


    """
        Lattice (L₁::_2VectorFloat, L₂::_2VectorFloat, G₁::_2VectorFloat, G₂::_2VectorFloat; originOffsetUV = [], originOffsetXY = [] )
    For internal use.
    Returns a lattice object taking in both real-space and reciprocal space vectors.
    If only one originOffset parameter has been set, that values is used to calculate the other one.  If neither originOffset parameter has been set, both are [0,0].
    """
    function Lattice(L₁::_2VectorFloat, L₂::_2VectorFloat, G₁::_2VectorFloat, G₂::_2VectorFloat; originOffsetUV = [], originOffsetXY = [], gridAlignment = CENTERALIGNMENT )

        # Calc coordinate transformation matrices
        transformationUVtoXY = [L₁[X] L₂[X]; L₁[Y] L₂[Y]]
        #Under certain symmetry, a singular exception occurs, so neec to offset just away from symmetry.
        if L₂ == [0,0]  # If a 1D lattice
            if L₁[X] == L₁[Y]
                transformationXYtoUV = inv( [L₁[X] eps(Float64); L₁[Y] 2*eps(Float64)] )
            else
                transformationXYtoUV = inv( [L₁[X] eps(Float64); L₁[Y] eps(Float64)] )
            end
        else
            transformationXYtoUV = inv(transformationUVtoXY)
        end

        # Calc originoffset
        if originOffsetUV != [] && originOffsetXY != []  # Both offsets have been set
            error("Can't set both a UV and XY offset in lattice.")
        elseif originOffsetUV == [] && originOffsetXY == []  # Neither offset has been set
            originOffsetUV = _2VectorFloat(0.0,0.0)
            originOffsetXY = _2VectorFloat(0.0,0.0)
        elseif originOffsetXY == []  # UV offset has been set
            originOffsetXY = transformationUVtoXY * originOffsetUV
        elseif originOffsetUV == []  # XY offset has been set
            originOffsetUV = transformationXYtoUV * originOffsetXY
        end


        return new(L₁,L₂,G₁,G₂, _2VectorFloat(originOffsetUV), _2VectorFloat(originOffsetXY), transformationUVtoXY, transformationXYtoUV, gridAlignment)
    end
end


# Define 2D lattice using two column vectors.
"""
    Lattice(L₁,L₂;originOffsetUV=[],originOffsetXY=[])
Returns a Lattice object with the periodicity of the L₁, L₂ lattice vectors.

If only one originOffset parameter has been set, that values is used to calculate the other one.  If neither originOffset parameter has been set, both are [0,0].
"""
function Lattice(L₁,L₂;originOffsetUV=[],originOffsetXY=[], gridAlignment = CENTERALIGNMENT)
    L₁ = convert(_2VectorFloat,L₁)
    L₂ = convert(_2VectorFloat,L₂)
    G₁, G₂ = calcReciprocals(L₁, L₂)
    return Lattice(L₁,L₂,G₁,G₂; originOffsetUV=originOffsetUV, originOffsetXY=originOffsetXY, gridAlignment=gridAlignment)
end

export calcReciprocals
"""
    calcReciprocals(L₁,L₂)
Calculates the reciprocal lattice vectors G₁, G₂ from the real-space lattice vectors, L₁, L₂.
"""
function calcReciprocals(L₁, L₂)
    @assert L₁ != [0,0] && L₂ != [0,0]
    G₁_G₂ = 2*π*transpose(inv([L₁ L₂]))
    return G₁_G₂[:,1], G₁_G₂[:,2]
end

"""
    calcReciprocals(L)
Calculates the reciprocal lattice vector, G₁, from the real-space lattice vector, L₁.  For a 1D lattice.
"""
function calcReciprocals(L)
    L = convert(_2VectorFloat,L)
    @assert L != _2VectorFloat(0,0)
    normL = norm(L)
    return convert(_2VectorFloat, 2*pi * L ./ norm(L)^2)
end

# If only given 1 vector, define lattice by setting one 2nd real and reciprocal space vectors to 0-vector
"""
    Lattice(L₁; originOffsetUV = [], originOffsetXY = [] )
Returns a Lattice object with the periodicity of the L₁ lattice vectors.  The L₂ lattice vector is [0,0].

If only one originOffset parameter has been set, that values is used to calculate the other one.  If neither originOffset parameter has been set, both are [0,0].
"""
function Lattice(L₁; originOffsetUV = [], originOffsetXY = [], gridAlignment = CENTERALIGNMENT)
    L₁ = convert(_2VectorFloat, L₁)
    G₁ = calcReciprocals( L₁ )
    L₂ = _2VectorFloat(0,0)
    G₂ = _2VectorFloat(0,0)
    return Lattice( L₁, L₂, G₁, G₂; originOffsetUV=originOffsetUV, originOffsetXY=originOffsetXY, gridAlignment )
end

# If only given a scalar, assume that is a 1D lattice only in x-direction
"""
    Lattice(L::Real; originOffsetUV = [], originOffsetXY = [] )
Returns a Lattice object with the equal to L. The L₁ lattice vector is [L,0]. The L₂ lattice vector is [0,0].

If only one originOffset parameter has been set, that values is used to calculate the other one.  If neither originOffset parameter has been set, both are [0,0].
"""
function Lattice(l₁::Real; originOffsetUV = [], originOffsetXY = [], gridAlignment = CENTERALIGNMENT)
    return Lattice( _2VectorFloat(l₁,0); originOffsetUV=originOffsetUV, originOffsetXY=originOffsetXY )
end


export is1D
"""
    is1D(lattice::Lattice)
Returns true if the lattice has only one periodicity.  False otherwise.
"""
function is1D(lattice::Lattice)
    return lattice.L₂ == [0,0]
end

export calcLatticeBoundaryLine
# Returns two parallel vectors for drawing the perimiter of the lattice unit cell.
"""
    calcLatticeBoundaryLine(lattice::Lattice)
Returns Xvector, Yvector
Returns 2 vectors representing the X and Y values of the unit cell lattice coordinates.
"""
function calcLatticeBoundaryLine(lattice::Lattice)
    L₁ = lattice.L₁
    L₂ = lattice.L₂
    Lx = [0, L₁[X], L₁[X]+L₂[X], L₂[X], 0] .+ lattice.originOffsetXY[X]
    Ly = [0, L₁[Y], L₁[Y]+L₂[Y], L₂[Y], 0] .+ lattice.originOffsetXY[Y]
    return Lx, Ly
end

export calcReciprocalLatticeBoundaryLine
# Returns two parallel vectors for drawing the perimiter of the reciprocal lattice unit cell.
"""
    calcReciprocalLatticeBoundaryLine(lattice::Lattice)
Returns Xvector, Yvector
Returns 2 vectors representing the X and Y values of the unit cell lattice coordinates.
"""
function calcReciprocalLatticeBoundaryLine(lattice::Lattice)
    G₁ = lattice.G₁
    G₂ = lattice.G₂
    Gx = [0, G₁[X], G₁[X]+G₂[X], G₂[X], 0]
    Gy = [0, G₁[Y], G₁[Y]+G₂[Y], G₂[Y], 0]
    return Gx, Gy
end


export getLatticeCenterCoordinates
"""
    getLatticeCenterCoordinates(lattice::Lattice)
Returns the X,Y coordinates of the center of the lattice
"""
function getLatticeCenterCoordinates(lattice::Lattice)
    return lattice.L₁*0.5 + lattice.L₂*0.5 + lattice.originOffsetXY
end

export convertUVtoXY
"""
    convertUVtoXY(lattice::Lattice, uv)
Converts the input U,V coordinates to equivalent X,Y coordinates
"""
function convertUVtoXY(lattice::Lattice, uv)
    return lattice.transformationUVtoXY * uv + lattice.originOffsetXY
end

export convertXYtoUV
"""
    convertXYtoUV(lattice::Lattice, uv)
Converts the input X,Y coordinates to equivalent U,V coordinates
"""
function convertXYtoUV(lattice::Lattice, xy)
    return lattice.transformationXYtoUV * (xy - lattice.originOffsetXY)
end
