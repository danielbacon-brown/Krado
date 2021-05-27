
mutable struct PlottingParameters
    # faces:
    color  # may be a vector or string
    alpha::Float64
    shade::Bool
    # edges:
    lineColor # may be a vector or string
    lineWidth::Float64
    lineStyle::String
    # vertex:
    marker::String

    function PlottingParameters(; color="gray", alpha=0.8, shade=true, lineColor="black", lineWidth=0.5, lineStyle="-", marker="," )
        return new(
            color,
            alpha,
            shade,
            lineColor,
            lineWidth,
            lineStyle,
            marker)
    end
end
