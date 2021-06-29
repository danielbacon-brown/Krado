export LayerStack

mutable struct LayerStack

    layers::Vector{<:LayerDefinition}

    function LayerStack(layers::Vector{<:LayerDefinition})
        return new(layers)
    end

end

function LayerStack()
    return LayerStack( Vector{<:LayerDefinition}() )
end

Base.first(layerStack::LayerStack) = first(layerStack.layers)
Base.last(layerStack::LayerStack) = last(layerStack.layers)
Base.length(layerStack::LayerStack) = length(layerStack.layers)
Base.getindex(layerStack::LayerStack, index::Integer) = getindex(layerStack.layers, index)
Base.iterate(layerStack::LayerStack) = iterate(layerStack.layers)
Base.iterate(layerStack::LayerStack, i::Integer) = iterate(layerStack.layers, i)
