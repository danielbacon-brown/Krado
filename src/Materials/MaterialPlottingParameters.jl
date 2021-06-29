export PlottingParameterCollection, addPlottingParameter!, getPlottingParameter

mutable struct PlottingParameterCollection
    dict::Dict{String,PlottingParameters}

    function PlottingParameterCollection(dict::Dict{String,PlottingParameters})
        return new(dict)
    end
end
# MaterialPlottingParameterCollection(dict::Dict{String,PlottingParameters}) = MaterialCollection(convert(Dict{String,AbstractMaterial}, dict) )
PlottingParameterCollection() = PlottingParameterCollection( Dict{String,PlottingParameters}() )

function addPlottingParameter!(plottingParameters::PlottingParameterCollection, name::String, params::PlottingParameters)
    plottingParameters.dict[name] = params
end

function getPlottingParameter(plottingParameters::PlottingParameterCollection, name::String)
    return plottingParameters.dict[name]
end

# Returns true if the dict contains the given name
function contains(plottingParameters::PlottingParameterCollection, name::String)
    return name in keys(plottingParameters.dict)
end

Base.getindex(plottingParamCol::PlottingParameterCollection, s::String) = getindex(plottingParamCol.dict, s)
Base.iterate(plottingParamCol::PlottingParameterCollection) = iterate(plottingParamCol.dict)
Base.iterate(plottingParamCol::PlottingParameterCollection, i::Integer) = iterate(plottingParamCol.dict, i)
