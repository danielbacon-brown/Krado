import Base.isless

# Connects a shape with the material that should be there
mutable struct Solid
    shape::Shape
    materialName::String
    priority::Float64
end

# Default to priority 1
function Solid( shape::Shape, materialName::String )
    return Solid( shape, materialName, 1)
end

function containsXY( solid::Solid, point)::Bool
    return containsXY( solid.shape, point )
end

# Establish sorting by priorty:
function isless(solid1::Solid, solid2::Solid)::Bool
    return solid1.priority < solid2.priority
end
