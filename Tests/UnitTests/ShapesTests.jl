@testset "Rectangle" begin
    # Test containsXY of Rectangle
    rect1 = Rectangle([0,0], [1,1])
    @test containsXY( rect1, [0, 0])
    @test containsXY( rect1, [0.45, 0])
    @test containsXY( rect1, [-0.45, 0])
    @test containsXY( rect1, [0, 0.45])
    @test containsXY( rect1, [0, -0.45])
    @test containsXY( rect1, [0.55, 0]) == false
    @test containsXY( rect1, [-0.55, 0]) == false
    @test containsXY( rect1, [0, 0.55]) == false
    @test containsXY( rect1, [0, -0.55]) == false
    rect2 = Rectangle([5,-4], [4,4])
    @test containsXY( rect2, [6.5,-4])
    @test containsXY( rect2, [3.5,-4])
    @test containsXY( rect2, [0,0]) == false
end;

@testset "Circle" begin
    # Test containsXY of circle
    circ1 = Circle([0,0], 1)
    @test containsXY( circ1, [0,0])
    @test containsXY( circ1, [0,0.9])
    @test containsXY( circ1, [0.9,0])
    @test containsXY( circ1, [0.9,0.9]) == false
    circ2 = Circle([10,2], 4)
    @test containsXY( circ2, [11,1])
    @test containsXY( circ2, [1,1]) == false
end;

@testset "Line segments" begin
    # Test segment intersection. Positive slope.
    vert1 = [0.0, 0.0]
    vert2 = [1.0, 2.0]
    @test rayIntersectsSegment(vert1, vert2, [1.5, 1.0]) == false
    @test rayIntersectsSegment(vert1, vert2, [0.5, 3.0]) == false
    @test rayIntersectsSegment(vert1, vert2, [0.5, -3.0]) == false
    @test rayIntersectsSegment(vert1, vert2, [-0.5, 1.0])
    @test rayIntersectsSegment(vert1, vert2, [0.5, 0.9]) == false
    @test rayIntersectsSegment(vert1, vert2, [0.5, 1.1])
    # Test segment intersection. Negative slope.
    vert1 = [0.0, 0.0]
    vert2 = [1.0, -2.0]
    @test rayIntersectsSegment(vert1, vert2, [1.5, -1.0]) == false
    @test rayIntersectsSegment(vert1, vert2, [0.5, -3.0]) == false
    @test rayIntersectsSegment(vert1, vert2, [0.5, 3.0]) == false
    @test rayIntersectsSegment(vert1, vert2, [-0.5, -1.0])
    @test rayIntersectsSegment(vert1, vert2, [0.5, -0.9]) == false
    @test rayIntersectsSegment(vert1, vert2, [0.5, -1.1])

    # Failed at polygon test level:
    vert1 = [2.0, 0.0]
    vert2 = [0.0, 0.2]
    testPoint = [2.0, 0.5]
    @test rayIntersectsSegment(vert1, vert2, testPoint) == false

    # Boolean operations
    bool1 = Rectangle([0,0],[3,3])
    bool2 = Circle([1,0],1.5)
    pt1 = [-1,0]
    pt2 = [0,0]
    pt3 = [2,0]
    pt4 = [0,5]

    union = UnionShape()
    append!(union, [bool1,bool2])
    @test containsXY(union, pt1) == true
    @test containsXY(union, pt2) == true
    @test containsXY(union, pt3) == true
    @test containsXY(union, pt4) == false

    intersection = IntersectionShape(bool1)
    push!(intersection,bool2)
    @test containsXY(intersection, pt1) == false
    @test containsXY(intersection, pt2) == true
    @test containsXY(intersection, pt3) == false
    @test containsXY(intersection, pt4) == false

    subtraction = SubtractionShape(bool1,[bool2])
    @test containsXY(subtraction, pt1) == true
    @test containsXY(subtraction, pt2) == false
    @test containsXY(subtraction, pt3) == false
    @test containsXY(subtraction, pt4) == false

    difference = DifferenceShape([bool1,bool2])
    @test containsXY(difference, pt1) == true
    @test containsXY(difference, pt2) == false
    @test containsXY(difference, pt3) == true
    @test containsXY(difference, pt4) == false

end;

@testset "Polygon" begin
    poly1 = Polygon([[0,0], [2,0], [0,2]])
    @test containsXY(poly1, [0.5,0.5])
    @test containsXY(poly1, [2.0,0.5]) == false
    @test containsXY(poly1, [-1,0.5]) == false
    @test containsXY(poly1, [3,0.5]) == false
    @test containsXY(poly1, [0.5,2]) == false

    poly2 = Polygon([[0,0], [1,-1], [0,1],[-1,-1]])
    @test containsXY(poly2, [0,0.5])
    @test containsXY(poly2, [0,-0.5]) == false
    @test containsXY(poly2, [0.5,-0.2])
    @test containsXY(poly2, [-0.5,-0.2])
    @test containsXY(poly2, [-0.5, 0.5]) == false

    # Testing equilateral triangle around 0
    Ly = 1.5 * cm
    triangle_length = 0.8 * Ly
    point1 = [triangle_length/2, -(triangle_length/2) / tand(30)]
    point2 = [-triangle_length/2, -(triangle_length/2) / tand(30)]
    point3 = [0, (triangle_length/2) / sind(60)]

    poly3 = Polygon( [point1, point2, point3])
    @test containsXY(poly3, [0,0.0])
    @test containsXY(poly3, [0.0001,0.0])
    @test containsXY(poly3, [0.0001,0.0001])
    @test containsXY(poly3, [0.014,0.0]) == false

    # TODO Test of polygon center:

end;
