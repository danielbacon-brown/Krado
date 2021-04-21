# Contains different kinds of coordinates (for plotting)
# RELATED FUNCTIONS:
# plot( SurfacePlotter, CoordinateData ) #
abstract type CoordinateData end


# Grid of coordinates where adjacent coordinates are adjacent in real space (for plotting)
mutable struct GridCoordinateData <: CoordinateData
	pointsX::Array{Float64,2}
	pointsY::Array{Float64,2}
	pointsZ::Array{Float64,2}
end

# Coordinates of triangulation of the surface.
mutable struct TriangulationCoordinateData <: CoordinateData

	# Coordinate of each vertex
	pointsX::Vector{Float64}
	pointsY::Vector{Float64}
	pointsZ::Vector{Float64}

	# n x 3 array of indices
	# Each row indicates a triangle
	triangleIndices::Array{Int64,2}
end

# Simple lists of coordinates
mutable struct PointCoordinateData <: CoordinateData
	pointsX::Vector{Float64}
	pointsY::Vector{Float64}
	pointsZ::Vector{Float64}
end

# Translate surface plot coordinates by parents position
function translate!( coordinateData::CoordinateData, position::Vector{Float64} )
	coordinateData.pointsX .+= position[XDIM]
	coordinateData.pointsY .+= position[YDIM]
	coordinateData.pointsZ .+= position[ZDIM]
end
