# MaterialCollection stores the pairing of materials with materials names.  Essentially just a dict, but easier to understand as it's own type.

mutable struct MaterialCollection
    dict::Dict{String,AbstractMaterial}

    function MaterialCollection(dict::Dict{String,AbstractMaterial})
        return new(dict)
    end
end
MaterialCollection(dict::Dict{String,T}) where T<:AbstractMaterial = MaterialCollection(convert(Dict{String,AbstractMaterial}, dict) )
MaterialCollection() = MaterialCollection( Dict{String,AbstractMaterial}() )

function addMaterial!(matCollection::MaterialCollection, name::String, mat::AbstractMaterial)
    matCollection.dict[name] = mat
end

function getMaterial(matCollection::MaterialCollection, name::String)
    return matCollection.dict[name]
end

# Returns true if the dict contains the given name
function contains(matCol::MaterialCollection, name::String)
    return name in keys(matCol.dict)
end
