# Creates a dictionary of: name => material import function
# function getFavoriteMaterialImporters(path::String)
# 
#     favorites = Dict{String,Function}()
# 
#     include(DIR*"\\"*path)
# 
#     # push!(favorites, "Ag"=> () -> importλnkTextMaterial("UserMaterials\\Ag_Johnson_and_Christy.txt"; scale=μm, skipRows=1, delimiter="\t"))
#     # push!(favorites, "Air"=> () -> importYAMLmaterial("MaterialDatabase\\data\\other\\mixed gases\\air\\Ciddor.yml"))
# 
#     return favorites
# end

# Import a single material from favorites
function importFavoriteMaterial!(matCol::MaterialCollection, materialImporters::Dict{<:Any,<:Any}, name::String)
    if name in keys(materialImporters)
        addMaterial!(matCol,name, materialImporters[name]() )
    end
end
# If there is no input MaterialCollection, create one.
function importFavoriteMaterial(materialImporters::Dict{<:Any,<:Any}, name::String)
    matCol = MaterialCollection()
    if name in keys(materialImporters)
        addMaterial!(matCol,name, materialImporters[name]() )
    end
    return matCol
end

# Import a series of materials from favorites
# names must be set-like and iterable
function importFavoriteMaterial!(matCol::MaterialCollection, materialImporters::Dict{<:Any,<:Any}, names)
    for name in names
        importFavoriteMaterial!(matCol, materialImporters, name)
    end
end
# If there is no material collection, create one
function importFavoriteMaterial(materialImporters::Dict{<:Any,<:Any}, names)
    matCol = MaterialCollection()
    importFavoriteMaterial!(matCol::MaterialCollection, materialImporters::Dict{<:Any,<:Any}, names)
    return matCol
end


# Returns a set of all material name used in the stack
function getMaterialsUsed(layerStack)::Set{String}
    nameSet = Set{String}()
    for layer in layerStack
        layerNames = getMaterialsUsed(layer)
        union!(nameSet, layerNames)
    end
    return nameSet
end

# Returns a set of all material names used in the layer
function getMaterialsUsed(layer::T) where T<:Union{UniformLayerDefinition, SemiInfiniteLayerDefinition}
    return Set{String}([layer.backgroundMaterialName])
end

# Returns a set of all material names used in the layer
function getMaterialsUsed(layer::PatternedLayerDefinition)::Set{String}
    nameSet = Set{String}()
    union!(nameSet, getMaterialsUsed(layer.layerPattern) )
    return nameSet
end

# Returns a set of all material names used in the layer
function getMaterialsUsed(pattern::LayerPattern)::Set{String}
    nameSet = Set{String}([pattern.backgroundMaterialName])
    for solid in pattern.solids
        push!(nameSet, solid.materialName )
    end
    return nameSet
end
