
# TODO: change r### and calcLatticeBoundaryLine, so that it is easier to understand which index corresponds to where in the lattice.  Maybe a dictionary.

function plot3VectorGrid(ax, grid::Array{_3VectorFloat,2}, params::PlottingParameters)
    x = map( r -> r[X], grid)
    y = map( r -> r[Y], grid)
    z = map( r -> r[Z], grid)
    ax.plot_surface( x, y, z, rstride=2,cstride=2, color=params.color,  alpha=params.alpha, shade=params.shade, edgecolor=params.lineColor, linewidth=params.lineWidth, linestyle=params.lineStyle )
end


function plot3Dsubstrate(ax, lattice::Lattice, layerStack::Vector{<:LayerDefinition}, materialParams::Dict{String, PlottingParameters}; scale=1)

    # lattice = simulationDefinition.lattice

    layer = first(layerStack)

    totalThickness = calcTotalThickness(layerStack)
    zLimits = getLayerStackPlotLimits(layerStack)

    zOuter = zLimits[1]
    zInner = 0.0

    # Follow lattice dimensions in X,Y
    latticeX, latticeY = calcLatticeBoundaryLine(lattice)

    # In r###, the first # is whether it has L₁, second is L₂, # is whether it is toward the middle of the stack or the outside of the stack.
    r111 = _3VectorFloat(latticeX[1], latticeY[1], zOuter)
    r211 = _3VectorFloat(latticeX[2], latticeY[2], zOuter)
    r221 = _3VectorFloat(latticeX[3], latticeY[3], zOuter)
    r121 = _3VectorFloat(latticeX[4], latticeY[4], zOuter)
    r112 = _3VectorFloat(latticeX[1], latticeY[1], zInner)
    r212 = _3VectorFloat(latticeX[2], latticeY[2], zInner)
    r222 = _3VectorFloat(latticeX[3], latticeY[3], zInner)
    r122 = _3VectorFloat(latticeX[4], latticeY[4], zInner)


    # Top left corner of matrix is lower X, higher Y
    top = _3VectorFloat[
            [r122] [r222];
            [r112] [r212] ] / scale
    south = _3VectorFloat[
            [r111] [r211];
            [r112] [r212] ] / scale
    north = _3VectorFloat[
            [r121] [r221];
            [r122] [r222] ] / scale
    west = _3VectorFloat[
            [r111] [r121];
            [r112] [r122] ] / scale
    east = _3VectorFloat[
            [r211] [r221];
            [r212] [r222] ] / scale

    plot3VectorGrid(ax, top, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, south, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, north, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, west, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, east, materialParams[layer.backgroundMaterialName])

end

function plot3Dsuperstrate(ax, lattice::Lattice, layerStack::Vector{<:LayerDefinition}, materialParams::Dict{String, PlottingParameters}; scale=1)

    # lattice = simulationDefinition.lattice

    layer = last(layerStack)

    totalThickness = calcTotalThickness(layerStack)
    zLimits = getLayerStackPlotLimits(layerStack)

    zOuter = zLimits[2]
    zInner = totalThickness

    # Follow lattice dimensions in X,Y
    latticeX, latticeY = calcLatticeBoundaryLine(lattice)

    # In r###, the first # is whether it has L₁, second is L₂, # is whether it is toward the middle of the stack or the outside of the stack.
    r111 = _3VectorFloat(latticeX[1], latticeY[1], zOuter)
    r211 = _3VectorFloat(latticeX[2], latticeY[2], zOuter)
    r221 = _3VectorFloat(latticeX[3], latticeY[3], zOuter)
    r121 = _3VectorFloat(latticeX[4], latticeY[4], zOuter)
    r112 = _3VectorFloat(latticeX[1], latticeY[1], zInner)
    r212 = _3VectorFloat(latticeX[2], latticeY[2], zInner)
    r222 = _3VectorFloat(latticeX[3], latticeY[3], zInner)
    r122 = _3VectorFloat(latticeX[4], latticeY[4], zInner)


    # Top left corner of matrix is lower X, higher Y
    bottom = _3VectorFloat[
            [r122] [r222];
            [r112] [r212] ] / scale
    south = _3VectorFloat[
            [r111] [r211];
            [r112] [r212] ] / scale
    north = _3VectorFloat[
            [r121] [r221];
            [r122] [r222] ] / scale
    west = _3VectorFloat[
            [r112] [r122];
            [r111] [r121] ] / scale
    east = _3VectorFloat[
            [r212] [r222];
            [r211] [r221] ] / scale

    plot3VectorGrid(ax, bottom, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, south, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, north, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, west, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, east, materialParams[layer.backgroundMaterialName])

end


function plot3Dlayer(ax, lattice::Lattice, layer::UniformLayerDefinition, zPosition::Real, materialParams::Dict{String, PlottingParameters}; scale=1)

    # lattice = simulationDefinition.lattice

    zLower = zPosition
    zUpper = (zPosition+layer.thickness)

    # Follow lattice dimensions in X,Y
    latticeX, latticeY = calcLatticeBoundaryLine(lattice)

    # In r###, the first # is whether it has L₁, second is L₂, # is whether it is toward the middle of the stack or the outside of the stack.
    r111 = _3VectorFloat(latticeX[1], latticeY[1], zLower)
    r211 = _3VectorFloat(latticeX[2], latticeY[2], zLower)
    r221 = _3VectorFloat(latticeX[3], latticeY[3], zLower)
    r121 = _3VectorFloat(latticeX[4], latticeY[4], zLower)
    r112 = _3VectorFloat(latticeX[1], latticeY[1], zUpper)
    r212 = _3VectorFloat(latticeX[2], latticeY[2], zUpper)
    r222 = _3VectorFloat(latticeX[3], latticeY[3], zUpper)
    r122 = _3VectorFloat(latticeX[4], latticeY[4], zUpper)


    # Top left corner of matrix is lower X, higher Y
    bottom = _3VectorFloat[
            [r121] [r221];
            [r111] [r211] ] / scale
    top = _3VectorFloat[
            [r122] [r222];
            [r112] [r212] ] / scale
    south = _3VectorFloat[
            [r111] [r211];
            [r112] [r212] ] / scale
    north = _3VectorFloat[
            [r121] [r221];
            [r122] [r222] ] / scale
    west = _3VectorFloat[
            [r112] [r122];
            [r111] [r121] ] / scale
    east = _3VectorFloat[
            [r212] [r222];
            [r211] [r221] ] / scale

    plot3VectorGrid(ax, bottom, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, top, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, south, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, north, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, west, materialParams[layer.backgroundMaterialName])
    plot3VectorGrid(ax, east, materialParams[layer.backgroundMaterialName])

end

function plot3Dlayer(ax, lattice::Lattice, layer::PatternedLayerDefinition, zPosition::Real, materialParams::Dict{String, PlottingParameters}; scale=1)
# function plotPatch3D(layer::PatternedLayerDefinition)

    # lattice = simulationDefinition.lattice

    zLower = zPosition
    zUpper = (zPosition+layer.thickness)

    # Follow lattice dimensions in X,Y
    latticeX, latticeY = calcLatticeBoundaryLine(lattice)

    # In r###, the first # is whether it has L₁, second is L₂, # is whether it is toward the middle of the stack or the outside of the stack.
    r111 = _3VectorFloat(latticeX[1], latticeY[1], zLower)
    r211 = _3VectorFloat(latticeX[2], latticeY[2], zLower)
    r221 = _3VectorFloat(latticeX[3], latticeY[3], zLower)
    r121 = _3VectorFloat(latticeX[4], latticeY[4], zLower)
    r112 = _3VectorFloat(latticeX[1], latticeY[1], zUpper)
    r212 = _3VectorFloat(latticeX[2], latticeY[2], zUpper)
    r222 = _3VectorFloat(latticeX[3], latticeY[3], zUpper)
    r122 = _3VectorFloat(latticeX[4], latticeY[4], zUpper)


    # Top left corner of matrix is lower X, higher Y
    bottom = _3VectorFloat[
            [r121] [r221];
            [r111] [r211] ] / scale
    top = _3VectorFloat[
            [r122] [r222];
            [r112] [r212] ] / scale
    south = _3VectorFloat[
            [r111] [r211];
            [r112] [r212] ] / scale
    north = _3VectorFloat[
            [r121] [r221];
            [r122] [r222] ] / scale
    west = _3VectorFloat[
            [r112] [r122];
            [r111] [r121] ] / scale
    east = _3VectorFloat[
            [r212] [r222];
            [r211] [r221] ] / scale

    plot3VectorGrid(ax, bottom, materialParams[layer.layerPattern.backgroundMaterialName])
    plot3VectorGrid(ax, top, materialParams[layer.layerPattern.backgroundMaterialName])
    plot3VectorGrid(ax, south, materialParams[layer.layerPattern.backgroundMaterialName])
    plot3VectorGrid(ax, north, materialParams[layer.layerPattern.backgroundMaterialName])
    plot3VectorGrid(ax, west, materialParams[layer.layerPattern.backgroundMaterialName])
    plot3VectorGrid(ax, east, materialParams[layer.layerPattern.backgroundMaterialName])


    # plot layerPattern:
    zRange = _2VectorFloat(zLower, zUpper)
    plot3DlayerPattern(ax, layer.layerPattern, zRange, materialParams; scale=scale)

end


function plot3DlayerPattern(ax, pattern::LayerPattern, zRange::TU2VectorReal, materialParams::Dict{String, PlottingParameters}; scale=1)

    for solid in pattern.solids
        plotPatch3D(ax, solid, zRange, materialParams; scale=scale)
    end
end

function plotPatch3D(ax, solid::Solid, zRange::TU2VectorReal, materialParams::Dict{String, PlottingParameters}; scale=1)

    patches = getPlotPatches3D(solid.shape, zRange; scale=scale)

    for patch in patches
        plot3VectorGrid(ax, patch, materialParams[solid.materialName])
    end
end

function getPlotPatches3D(rect::Rectangle, zRange::TU2VectorReal; scale=1)

    patches = Vector{Array{_3VectorFloat,2}}(undef, 6)

    # In r###, the first # is whether it has L₁, second is L₂, # is whether it is toward the middle of the stack or the outside of the stack.
    x = [rect.center[X]-rect.lengths[X]/2, rect.center[X]+rect.lengths[X]/2]
    y = [rect.center[Y]-rect.lengths[Y]/2, rect.center[Y]+rect.lengths[Y]/2]
    r111 = _3VectorFloat(x[1], y[1], zRange[1])
    r211 = _3VectorFloat(x[2], y[1], zRange[1])
    r221 = _3VectorFloat(x[2], y[2], zRange[1])
    r121 = _3VectorFloat(x[1], y[2], zRange[1])
    r112 = _3VectorFloat(x[1], y[1], zRange[2])
    r212 = _3VectorFloat(x[2], y[1], zRange[2])
    r222 = _3VectorFloat(x[2], y[2], zRange[2])
    r122 = _3VectorFloat(x[1], y[2], zRange[2])


    # Top left corner of matrix is lower X, higher Y
    patches[1] = _3VectorFloat[
            [r121] [r221];
            [r111] [r211] ] / scale
    patches[2] = _3VectorFloat[
            [r122] [r222];
            [r112] [r212] ] / scale
    patches[3] = _3VectorFloat[
            [r111] [r211];
            [r112] [r212] ] / scale
    patches[4] = _3VectorFloat[
            [r121] [r221];
            [r122] [r222] ] / scale
    patches[5] = _3VectorFloat[
            [r112] [r122];
            [r111] [r121] ] / scale
    patches[6] = _3VectorFloat[
            [r212] [r222];
            [r211] [r221] ] / scale

    return patches
end

function getPlotPatches3D(circ::Circle, zRange::TU2VectorReal; scale=1)

    numPatches = 40

    # Polar angle, θ
    # There is one more angle than there are patches
	θᵢ = Vector(LinRange(0.0,2pi,numPatches+1))
    xᵢ = cos.(θᵢ)*circ.radius .+ circ.center[X]
    yᵢ = sin.(θᵢ)*circ.radius .+ circ.center[Y]

    patches = Vector{Array{_3VectorFloat,2}}(undef, numPatches)

    for iPatch in 1:numPatches
        r11 = _3VectorFloat(xᵢ[iPatch], yᵢ[iPatch], zRange[1])
        r12 = _3VectorFloat(xᵢ[iPatch], yᵢ[iPatch], zRange[2])
        r21 = _3VectorFloat(xᵢ[iPatch+1], yᵢ[iPatch+1], zRange[1])
        r22 = _3VectorFloat(xᵢ[iPatch+1], yᵢ[iPatch+1], zRange[2])
        patches[iPatch] = _3VectorFloat[ [r11] [r21]; [r12] [r22] ]/scale
    end

    return patches
end


function getPlotPatches3D(poly::Polygon, zRange::TU2VectorReal; scale=1)

    numPatches = length(poly.vertices)

    xᵢ = [poly.vertices[i][X] for i in 1:numPatches] .+ poly.offset[X]
    yᵢ = [poly.vertices[i][Y] for i in 1:numPatches] .+ poly.offset[Y]

    patches = Vector{Array{_3VectorFloat,2}}(undef, numPatches)

    for iPatch in 1:(numPatches-1)
        # @show iPatch
        r11 = _3VectorFloat(xᵢ[iPatch], yᵢ[iPatch], zRange[1])
        r12 = _3VectorFloat(xᵢ[iPatch], yᵢ[iPatch], zRange[2])
        r21 = _3VectorFloat(xᵢ[iPatch+1], yᵢ[iPatch+1], zRange[1])
        r22 = _3VectorFloat(xᵢ[iPatch+1], yᵢ[iPatch+1], zRange[2])
        patches[iPatch] = _3VectorFloat[ [r11] [r21]; [r12] [r22] ]/scale
    end
    r11 = _3VectorFloat(xᵢ[numPatches], yᵢ[numPatches], zRange[1])
    r12 = _3VectorFloat(xᵢ[numPatches], yᵢ[numPatches], zRange[2])
    r21 = _3VectorFloat(xᵢ[1], yᵢ[1], zRange[1])
    r22 = _3VectorFloat(xᵢ[1], yᵢ[1], zRange[2])
    patches[end] = _3VectorFloat[ [r11] [r21]; [r12] [r22] ]/scale

    return patches
end


function getPlotPatches3D(parentShape::T, zRange::TU2VectorReal; scale=1) where T<:Union{UnionShape, IntersectionShape, DifferenceShape}

    patches = Array{_3VectorFloat,2}[]
    for shape in parentShape.shapes
        append!(patches, getPlotPatches3D(shape, zRange; scale=scale) )
    end
    return patches
end

function getPlotPatches3D(parentShape::SubtractionShape, zRange::TU2VectorReal; scale=1)

    patches = Array{_3VectorFloat,2}[]
    append!(patches, getPlotPatches3D(parentShape.baseShape, zRange; scale=scale) )
    for shape in parentShape.shapes
        append!(patches, getPlotPatches3D(shape, zRange; scale=scale) )
    end
    return patches
end
