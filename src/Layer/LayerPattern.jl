# This calculates the permittivity and permeability across space (x,y)
# These can be nested, forming a "pipeline" for calculating (x,y)
# For example, can do a boolean union/intersection of two solids

abstract type AbstractLayerPattern end

# Calculates the permittivity for a set of solids according to the priorty of each solid.
mutable struct LayerPattern <: AbstractLayerPattern
    solids::Vector{Solid}
    backgroundMaterialName::String

    function LayerPattern( solids::Vector{Solid}, backgroundMaterialName::String)
        solids = sort(solids, by = x -> x.priority)
        return new(solids, backgroundMaterialName)
    end
end

function LayerPattern( backgroundMaterialName::String)
    return LayerPattern(Solid[], backgroundMaterialName)
end
function LayerPattern( solid::Solid, backgroundMaterialName::String)
    return LayerPattern([solid], backgroundMaterialName)
end

function addSolid(spatcalc::LayerPattern, solid::Solid)
    push!(spatcalc.solids, solid)
    sort!(spatcalc)
    return
end

function getMaterialAtPosition( spat::LayerPattern, point::_2VectorFloat )::String
    # Sort solids based on priority if needed, then calculate.
    if !issorted(spat)
        sort!(spat)
    end

    # Use first material than matches and the background material of there's no match
    for solid in spat.solids
        if containsXY(solid, point)
            return solid.materialName
        end
    end
    return spat.backgroundMaterialName
end
getMaterialAtPosition( spat::LayerPattern, point ) = getMaterialAtPosition( spat, _2VectorFloat(point) )

function getMaterialAtPosition( spat::LayerPattern, positionGrid::PositionGridXY{N} )::Array{String,N} where {N}

    return map(position -> getMaterialAtPosition(spat,position), positionGrid.positions)
end


import Base.sort!
import Base.issorted

function sort!(spat::LayerPattern)
    return sort!(spat.solids)
end

function issorted(spat::LayerPattern)
    return issorted(spat.solids)
end

#Returns the ϵμ corresponding to input position
function getϵμAtPosition( spat::LayerPattern, position::_2VectorFloat, materials::MaterialCollection, wavenumber::Wavenumber)
    matName = getMaterialAtPosition(spat, position )
    return calc_ϵμ( getMaterial(materials,matName), wavenumber, position)
end
getϵμAtPosition( spat::LayerPattern, position, materials::MaterialCollection, wavenumber::Wavenumber) = getϵμAtPosition( spat, _2VectorFloat(position), materials::MaterialCollection, wavenumber::Wavenumber)

# Declare methods for just getting a single characteristic (ϵ or μ):
getϵAtPosition( spat::LayerPattern, position::_2VectorFloat, materials::MaterialCollection, wavenumber::Wavenumber) = getϵμAtPosition( spat, position, materials, wavenumber)[1]
getϵAtPosition( spat::LayerPattern, position, materials::MaterialCollection, wavenumber::Wavenumber) = getϵμAtPosition( spat, _2VectorFloat(position), materials, wavenumber)[1]
getμAtPosition( spat::LayerPattern, position::_2VectorFloat, materials::MaterialCollection, wavenumber::Wavenumber) = getϵμAtPosition( spat, position, materials, wavenumber)[2]
getμAtPosition( spat::LayerPattern, position, materials::MaterialCollection, wavenumber::Wavenumber) = getϵμAtPosition( spat, _2VectorFloat(position), materials, wavenumber)[2]


# For array of positions
function getϵμAtPosition( spat::LayerPattern, positionGrid::PositionGridXY, materials::MaterialCollection, wavenumber::Wavenumber) where {N}
    return map( position -> getϵμAtPosition( spat, position, materials, wavenumber), positionGrid.positions )
end
function getϵAtPosition( spat::LayerPattern, positionGrid::PositionGridXY, materials::MaterialCollection, wavenumber::Wavenumber) where {N}
    return map( position -> getϵAtPosition( spat, position, materials, wavenumber), positionGrid.positions )
end
function getμAtPosition( spat::LayerPattern, positionGrid::PositionGridXY, materials::MaterialCollection, wavenumber::Wavenumber) where {N}
    return map( position -> getμAtPosition( spat, position, materials, wavenumber), positionGrid.positions )
end
