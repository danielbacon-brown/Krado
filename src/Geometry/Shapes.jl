# include("RayCastings.jl")
# import RayCastings

# Shape just includes raw geometry data.  Enough to see if something is inside or not.  Shape can contain other shapes (array or rotation matrix)
# containsXY(shape, coordinates)
abstract type Shape end

mutable struct Rectangle <: Shape
    center::_2VectorFloat
    lengths::_2VectorFloat

    function Rectangle(center::_2VectorFloat, lengths::_2VectorFloat)
        new(center, lengths)
    end
end
Rectangle( center, lengths) = Rectangle( _2VectorFloat(center), _2VectorFloat(lengths))

function containsXY(rect::Rectangle, testPoint)::Bool
    return all( broadcast( abs, testPoint-rect.center) .< 0.5*rect.lengths )
end


mutable struct Circle <: Shape
    center::_2VectorFloat
    radius::Float64

    function Circle(center::_2VectorFloat, radius::Float64)
        @assert radius > 0
        return new(center, radius)
    end
end
Circle( center, radius) = Circle(_2VectorFloat(center), Float64(radius) )


function containsXY(circ::Circle, testPoint)::Bool
    return sum( (testPoint-circ.center).^2 ) < circ.radius^2
end


mutable struct Polygon <: Shape
    # Coordinates is a list of x,y pairs
    center::_2VectorFloat
    vertices::Vector{_2VectorFloat}

    function Polygon( center::_2VectorFloat, vertices::Vector{_2VectorFloat} )
        @assert length(vertices) > 2
        return new(center, vertices)
    end
end
Polygon(center, vertices) = Polygon(_2VectorFloat(center), convert(Vector{_2VectorFloat}, vertices) )

# TODO: switch to using default params
# Default [0,0] offset
Polygon(vertices) = Polygon(_2VectorFloat(0,0),vertices)


# Draws a ray to the right of test point.  If this ray intersects the line segment of the two vertices, return True.  Otherwise False.
# Adapted from https://rosettacode.org/wiki/Ray-casting_algorithm#Julia
function rayIntersectsSegment(vert1, vert2, test)::Bool

    #Ensure first is below second
    if vert1[Y] > vert2[Y]
        vertLow, vertHigh = vert2, vert1
    else
        vertLow, vertHigh = vert1, vert2
    end

    # Ensure point is not exactly same Y-value as endpoints
    if (test[Y] == vertLow[Y]) || (test[Y] == vertHigh[Y])
        test = _2VectorFloat( test[X], nextfloat(test[Y]) )
    end
    # Ensure point is not exactly same X-value as endpoints
    if (test[X] == vertLow[X]) || (test[X] == vertHigh[X])
        test = _2VectorFloat( nextfloat(test[X]), test[Y] )
    end

    # No intersection if the segment is entirely above or below segment or to the left
    if (test[Y] < vertLow[Y]) || (test[Y] > vertHigh[Y]) || (test[X] > max(vertLow[X],vertHigh[X]))
        return false
    end
    # Guaranteed intersection if the segment is otherwise completely to the right
    if (test[X] < vertLow[X]) && (test[X] < vertHigh[X])
        return true
    end

    # If vertical line (and thus colinear with test point), say no intersection
    if vertLow[X] == vertHigh[X]
        return false
    end

    slopeVert = (vertHigh[Y] - vertLow[Y])/(vertHigh[X] - vertLow[X])
    slopeTest = (test[Y] - vertLow[Y] )/(test[X] - vertLow[X])
    return slopeTest > slopeVert

end



# Draws a ray in +x direction.  if the number of intersections is odd, it is inside the polygon
function containsXY(poly::Polygon, testPoint)::Bool
    testPoint = testPoint - poly.center

    numIntersections = 0
    # Test vertex pairs
    iₚ = 1
    for iₚ = 1:(length(poly.vertices)-1)
        if rayIntersectsSegment( poly.vertices[iₚ], poly.vertices[iₚ+1], testPoint)
            numIntersections += 1
        end
    end
    # Same test for first and last elements
    if rayIntersectsSegment( poly.vertices[1], poly.vertices[length(poly.vertices)], testPoint)

        numIntersections += 1
    end
    return isodd(numIntersections)
end



# Boolean functions of shapes.  A collection of shapes that functions as a shape itself
# All instances have shapes::Vector{Shape}
abstract type BooleanShape <: Shape end

import Base.push!
import Base.append!
# Adding shapes, individually or as vector
function push!(booleanShape::BooleanShape, shape::Shape)
    push!(booleanShape.shapes, shape)
end
function append!(booleanShape::BooleanShape, shapes)
    append!(booleanShape.shapes, shapes)
end


# The union of a set of shapes, that functions as a shape itself
mutable struct UnionShape <: BooleanShape
    shapes::Vector{Shape}
    function UnionShape(shapes::Vector{T}) where {T<:Shape}
        return new(shapes)
    end
end

# Can create with an empty vector, or with a single element
function UnionShape()
    return UnionShape(Shape[])
end
function UnionShape(shape::Shape)
    return UnionShape([shape])
end

function containsXY(unionShape::UnionShape, testPoint)::Bool
    return any([containsXY(shape, testPoint) for shape in unionShape.shapes])
end
# containsXY(unionShape::UnionShape, testPoint) = containsXY(unionShape, _2VectorFloat(testPoint))


# The intersection of a set of shapes, that functions as a shape itself.
# A point is in this shape if it is in all of the children shapes.
mutable struct IntersectionShape <: BooleanShape
    shapes::Vector{Shape}
    function IntersectionShape(shapes::Vector{T}) where {T<:Shape}
        return new(shapes)
    end
end

# Can create with an empty vector, or with a single element
function IntersectionShape()
    return IntersectionShape(Shape[])
end
function IntersectionShape(shape::Shape)
    return IntersectionShape([shape])
end

function containsXY(intersectionShape::IntersectionShape, testPoint::_2VectorFloat)::Bool
    return all([containsXY(shape, testPoint) for shape in intersectionShape.shapes])
end
containsXY(intersectionShape::IntersectionShape, testPoint) = containsXY(intersectionShape, _2VectorFloat(testPoint) )


# The difference of a set of shapes, that functions as a shape itself
# Difference defined as: base shape minus all other shapes.
mutable struct SubtractionShape <: BooleanShape
    baseShape::Shape
    shapes::Vector{Shape}
    function SubtractionShape(baseShape::Shape, shapes::Vector{T}) where {T<:Shape}
        return new(baseShape, shapes)
    end
end

# Can create with an empty vector, or with a single element
function SubtractionShape(baseShape::Shape)
    return SubtractionShape(baseShape,Shape[])
end
function SubtractionShape(baseShape::Shape, shape::Shape)
    return SubtractionShape(baseShape,[shape])
end

# Returns true if base shape contains point, but no shape in list does
function containsXY(subtractionShape::SubtractionShape, testPoint)::Bool
    inBase::Bool = containsXY(subtractionShape.baseShape, testPoint)
    inRest::Bool = any([containsXY(shape, testPoint) for shape in subtractionShape.shapes])
    return  inBase && !inRest
end


# The intersection of a set of shapes, that functions as a shape itself.
# A point is inside this shape if it is inside exactly one of the children shapes.
mutable struct DifferenceShape <: BooleanShape
    shapes::Vector{Shape}
    function DifferenceShape(shapes::Vector{T}) where {T<:Shape}
        return new(shapes)
    end
end

# Can create with an empty vector, or with a single element
function DifferenceShape()
    return IntersectionShape(Shape[])
end
function DifferenceShape(shape::Shape)
    return IntersectionShape([shape])
end

function containsXY(differenceShape::DifferenceShape, testPoint)::Bool
    return sum([containsXY(shape, testPoint) for shape in differenceShape.shapes]) == true
end
