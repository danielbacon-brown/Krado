# Defines the layer according to a prebuilt ϵ and μ grid
# Not the preferred way of defining a pattern.  User should be cautious and ensure that this is consistent with the lattice
mutable struct GridLayerPattern <: AbstractLayerPattern
    materialNameGrid::Array{String, 2}

    function GridLayerPattern(materialNameGrid::Array{String, 2})
        return new(materialNameGrid)
    end
end

# Returns ϵ,μ grids
# Note: this method returns the same result, regardless of the input positions.
function getϵμAtPosition( pattern::GridLayerPattern, positions::PositionGridXY, materials::MaterialCollection, wavenumber::Wavenumber) where {N}
    return map( materialName -> calc_ϵμ( getMaterial(materials,materialName), wavenumber), pattern.materialNameGrid )
end
function getϵAtPosition( pattern::GridLayerPattern, positions::PositionGridXY, materials::MaterialCollection, wavenumber::Wavenumber) where {N}
    return map( materialName -> calc_ϵ( getMaterial(materials,materialName), wavenumber), pattern.materialNameGrid )
end
function getμAtPosition( pattern::GridLayerPattern, positions::PositionGridXY, materials::MaterialCollection, wavenumber::Wavenumber) where {N}
    return map( materialName -> calc_μ( getMaterial(materials,materialName), wavenumber), pattern.materialNameGrid )
end
