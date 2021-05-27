mutable struct LayerStack

    layers::Vector{<:LayerDefinition}

    function LayerStack(layers::Vector{<:LayerDefinition})
        return new(layers)
    end

end

function LayerStack()
    return LayerStack( Vector{<:LayerDefinition}() )
end

Base.first(layerStack::LayerStack) = Base.first(layerStack.layers)
Base.last(layerStack::LayerStack) = Base.last(layerStack.layers)
Base.length(layerStack::LayerStack) = Base.length(layerStack.layers)
Base.getindex(layerStack::LayerStack, index::Integer) = Base.getindex(layerStack.layers, index)
