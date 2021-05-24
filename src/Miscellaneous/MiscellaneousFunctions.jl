

# Setting so that new plots get plotted in a standalone window.
pygui(true)


# Used to determine how the points are positioned in the unit cell.
abstract type GridAlignment end
mutable struct LeftAlignment <: GridAlignment end
mutable struct CenterAlignment <: GridAlignment end
mutable struct RightAlignment <: GridAlignment end
const global LEFTALIGNMENT = LeftAlignment()
const global CENTERALIGNMENT = CenterAlignment()
const global RIGHTALIGNMENT = RightAlignment()


# calculate the Fresnel coefficients between interfaces
function calcFresnelCoefficients(n₁, n₂, θ₁)
    θᵢ = θ₁
    θₜ = asin(n₁/n₂ * sin(θ₁))
    rₛ = ( n₁*cos(θᵢ) - n₂*cos(θₜ) )/( n₁*cos(θᵢ) + n₂*cos(θₜ) )
    rₚ = ( n₂*cos(θᵢ) - n₁*cos(θₜ) )/( n₂*cos(θᵢ) + n₁*cos(θₜ) )
    tₛ = 2*n₁*cos(θᵢ) / (n₁*cos(θᵢ) + n₂*cos(θₜ))
    tₚ = 2*n₁*cos(θᵢ) / (n₂*cos(θᵢ) + n₁*cos(θₜ))
    return tₛ, tₚ, rₛ, rₚ
end

# Generates a vector of evenly spaced values between 0 and 1.  Used to uniformly sample a periodic function.  "LeftAlignment" starts at zero
function unitLinspaceLeft(len::Integer)::Vector{Float64}
    @assert len > 0
    return convert(Vector{Float64}, range(0, 1, length=len+1)[1:end-1] )
end
function unitLinspace(gridAlignment::LeftAlignment, len::Integer)::Vector{Float64}
    return unitLinspaceLeft(len)
end
function unitLinspace(gridAlignment::CenterAlignment, len::Integer)::Vector{Float64}
    return unitLinspaceLeft(len) .+ 0.5/len
end
function unitLinspace(gridAlignment::RightAlignment, len::Integer)::Vector{Float64}
    return unitLinspaceLeft(len) .+ 1/len
end


# Creates a grid of tuples for every position on a grid with the given shape.
# INPUT: (int X, int Y) or _2VectorInt(X,Y)
# OUTPUT: [ (X₁,Y₁) (X₂,Y₁) ... (Xₙ,Y₁) ;
#           (X₁,Y₂) (X₂,Y₂) ... (Xₙ,Y₂) ;
#           ...
#           (X₁,Yₙ) (X₂,Yₙ) ... (Xₙ,Yₙ) ]
function getGridIndices(gridShape)::Array{Tuple{Int64,Int64},2}
    indices = collect( Base.product(1:gridShape[1], 1:gridShape[2]) )
    return reshape(indices, Tuple(gridShape) )
end

function unitize(𝐯)
    return 𝐯/norm(𝐯)
end



# Returns a 2D vector that has the same direction as L, but with an inverse length
function vectorInverse(v::_2VectorFloat)
    return unitize(v) * 1/norm(v)
end
vectorInverse(v)=unitize(_2VectorFloat(v))


# Safely returns half the integer
function half(n::T)::T where T<:Integer
    @assert (n % 2) == 0
    return round(T, n/2)
end
# Safely returns quarter the integer
function quarter(n::T)::T where T<:Integer
    @assert (n % 4) == 0
    return round(T, n/4)
end

# Returns +1 if true, -1 if negative
function bool2posNeg(isPositive)
    return isPositive ? (+1) : (-1)
end

# Returns the X,Y components of the given vector
function getXY(v::_3VectorFloat)::_2VectorFloat
    return _2VectorFloat(v[X:Y])
end
# getXY(v::TU3VectorReal) = getXY(_3VectorFloat(v))
function getXY(v::_3VectorComplex)::_2VectorComplex
    return _2VectorComplex(v[X:Y])
end
# getXY(v::TU3VectorReal) = getXY(_3VectorComplex(v))
function getXY(v::_3VectorInt)::_2VectorInt
    return _2VectorInt(v[X:Y])
end
getXY(v) = v[X:Y]

# Convert vacuum wavelength to vacuum wavenumber
function λ₀2k₀(λ₀)
    return 2*π/λ₀
end
# Convert vaccum wavenumber to vacuum wavelength
function k₀2λ₀(k₀)
    return 2*π/k₀
end


# k must be 3-vector.  A must be 2-vector
function fieldSPtoFieldXYZ(k, A)
    @assert length(k) == 3
    @assert length(A) == 2
    ŝ,p̂ = calcŝp̂(k)
    return ŝ*A[S] + p̂*A[P]
end


function fieldXYZtoFieldSP(k, E)
    @assert length(k) == 3
    @assert length(E) == 3
    ŝ,p̂ = calcŝp̂(k)
    return _2VectorComplex(ŝ ⋅ E, p̂ ⋅ E)
end

# TODO: combine this with above function
function Exyz2Esp(k::_3VectorFloat, E::_3VectorComplex)
    ŝ,p̂ = calcŝp̂(k)
    return _2VectorComplex(ŝ ⋅ E, p̂ ⋅ E)
end
Exyz2Esp(k::TU3VectorReal, E::TU3VectorComplex) = Exyz2Esp(_3VectorFloat(k), _3VectorComplex(E))

# Creates a float 3-vector corresponding to wavelength perpendicular to the surface.  Usually then rotated using rotateVector functions
function normalIncidenceKvector(wavenumber::Wavenumber)::_3VectorFloat
    return _3VectorFloat(0,0,getk₀(wavenumber))
end

# Rotates the incident 3-vector according to the azimuthal, θ, and zenith angle, ϕ.  Radians.
# Rotate ϕ around y-axis, then θ around z-axis
# rotY = [cosd(ϕ)     0       sind(ϕ);
#         0           1       0;
#         -sind(ϕ)    0       cosd(ϕ)]
# rotZ = [cosd(θ)     -sind(θ)    0;
#         sind(θ)     cosd(θ)     0;
#         0           0           1]
# returns rotZ * rotY * k
function rotateVector(k, θ::Real, ϕ::Real)::_3VectorFloat
    k = convert(_3VectorFloat,k)
    rotY = [cos(ϕ)     0       sin(ϕ);
            0           1       0;
            -sin(ϕ)    0       cos(ϕ)]
    rotZ = [cos(θ)     -sin(θ)    0;
            sin(θ)     cos(θ)     0;
            0           0           1]
    return rotZ * (rotY * k)
end
rotateVectorDegrees(k, θ::Real, ϕ::Real) = rotateVector(k, deg2rad(θ), deg2rad(ϕ))



function cartesian2SphericalCoordinates(kXYZ) ::Tuple{Float64, Float64, Float64}
    ρ = norm(kXYZ)
    θ = atan(Float64(kXYZ[Y]), Float64(kXYZ[X]) )
    ϕ = acos( kXYZ[Z]/ ρ)
    return ρ, θ, ϕ
end

# Uses radians.
function spherical2CartesianCoordinates(ρ::Real, θ::Real, ϕ::Real)
    x = ρ*sin(ϕ)*cos(θ)
    y = ρ*sin(ϕ)*sin(θ)
    z = ρ*cos(ϕ)
    return _3VectorFloat(x, y, z)
end
# Uses degrees.
spherical2CartesianCoordinatesDegrees(ρ::Real, θ::Real, ϕ::Real) = spherical2CartesianCoordinates( deg2rad(ρ), deg2rad(θ) )



# Calculate vector for s-polarization for the given k-vector
# For a k-vector only in z-direction, then the s-vector is parallel to +x̂
# ŝ = n̂×k / |n̂×k|
function calcŝ(k::_3VectorFloat)::_3VectorFloat
    crossprod = k × ẑ
    if crossprod == _3VectorFloat(0,0,0)
        return _3VectorFloat(1,0,0)
    end
    return unitize(crossprod)
end
# If it takes a generic, convert it to the specific type first
calcŝ(k::TU3VectorReal) = calcŝ(_3VectorFloat(k))

function calcŝ(k::_3VectorComplex)::_3VectorComplex
    crossprod = k × ẑ
    if crossprod == _3VectorComplex(0,0,0)
        return _3VectorComplex(1,0,0)
    end
    return unitize(crossprod)
end
# If it takes a generic, convert it to the specific type first
calcŝ(k::TU3VectorComplex) = calcŝ(_3VectorComplex(k))


# Calculate vector for p-polarization for the given k-vector
# p̂ = k×ŝ / |k×p̂|
function calcp̂(k::_3VectorFloat)::_3VectorFloat
    ŝ = calcŝ(k)
    return unitize(k × ŝ)
end
calcp̂(k::TU3VectorReal) = calcp̂(_3VectorFloat(k))

function calcp̂(k::_3VectorFloat, ŝ::_3VectorFloat)::_3VectorFloat
    return unitize(k × ŝ)
end
calcp̂(k::TU3VectorReal,ŝ::TU3VectorReal)=calcp̂(_3VectorFloat(k),_3VectorFloat(ŝ))


function calcp̂(k::_3VectorComplex)::_3VectorComplex
    ŝ = calcŝ(k)
    return unitize(k × ŝ)
end
calcp̂(k::TU3VectorComplex) = calcp̂(_3VectorComplex(k))

function calcp̂(k::_3VectorComplex, ŝ::_3VectorComplex)::_3VectorComplex
    return unitize(k × ŝ)
end
calcp̂(k::TU3VectorComplex,ŝ::TU3VectorComplex)=calcp̂(_3VectorComplex(k),_3VectorComplex(ŝ))


# Returns unit-vector of TE (s) and TM (p) polarization
# ŝ = n̂×k / |n̂×k|
# p̂ = k×ŝ / |k×p̂|
function calcŝp̂(k::_3VectorFloat)::Tuple{_3VectorFloat,_3VectorFloat}
    ŝ = calcŝ(k)
    return ŝ, calcp̂(k, ŝ)
end
calcŝp̂(k::TU3VectorReal) = calcŝp̂(_3VectorFloat(k))


function calcŝp̂(k::_3VectorComplex)::Tuple{_3VectorComplex,_3VectorComplex}
    ŝ = calcŝ(k)
    return ŝ, calcp̂(k, ŝ)
end
calcŝp̂(k::TU3VectorComplex) = calcŝp̂(_3VectorComplex(k))


# Takes a grid of positions and returns a vector of x-coordinates and a paired vector of y-coordinates. In the vectors, X changes more frequently.
function linearizePositionGrid(posGrid::Array{_2VectorFloat,2})::Tuple{ Vector{Float64}, Vector{Float64}}

    gridSize = size(posGrid)
    numElements = prod(gridSize)
    xCoords = Vector{Float64}(undef, numElements)
    yCoords = Vector{Float64}(undef, numElements)

    for i_x = 1:gridSize[X]
        for i_y = 1:gridSize[Y]
            xCoords[i_x + (i_y-1)*gridSize[X]] = posGrid[i_x,i_y][X]
            yCoords[i_x + (i_y-1)*gridSize[X]] = posGrid[i_x,i_y][Y]
        end
    end
    return xCoords, yCoords
end

# Takes a grid of strings referring to materialNames and reshapes to a vector. Designed to match linearizePositionGrid().
function linearizeMaterialNameGrid(materialNameGrid::Array{String,2})::Vector{String}

    gridSize = size(materialNameGrid)
    numElements = prod(gridSize)
    materialNames = Vector{String}(undef, numElements)


    for i_x = 1:gridSize[X]
        for i_y = 1:gridSize[Y]
            materialNames[i_x + (i_y-1)*gridSize[X]] = materialNameGrid[i_x,i_y]
        end
    end
    return materialNames
end

# Takes a grid of positions and returns a grid of x-coord and y-coord
function extractPositionGridComponents(posGrid::Array{_2VectorFloat,2}) ::Tuple{Array{Float64,2},Array{Float64,2}}
    gridSize = size(posGrid)
    xGrid = Array{Float64,2}(undef, gridSize)
    yGrid = Array{Float64,2}(undef, gridSize)

    for i_x = 1:gridSize[X]
        for i_y = 1:gridSize[Y]
            xGrid[i_x, i_y] = posGrid[i_x, i_y][X]
            yGrid[i_x, i_y] = posGrid[i_x, i_y][Y]
        end
    end

    return xGrid, yGrid
end




# Convert n, k into ϵ.  # Assumes that μ==1
function convert_n2ϵ(n::T)::T where T<:Number
    nReal, nImag = real(n), imag(n)
    ϵReal = nReal^2 - nImag^2
    ϵImag = 2*nReal*nImag
    return ϵReal + 1im*ϵImag
end


# Convert ϵ', ϵ'' into n,k.  # Assumes that μ==1
# k is negative when ϵ'' is negative.
function convert_ϵ2n(ϵ::T)::T where T<:Number
    return sqrt(ϵ)
end


# Convert ϵ', ϵ'' into n,k.  #  μ is not necessarily 1
# k is negative when ϵ'' is negative.
function convert_ϵμ2n(ϵ::T1, μ::T2) where {T1<:Number, T2<:Number}
    return sqrt(ϵ*μ)
end


# Converts relative ϵ μ to relative impedance z
function convert_ϵμ2z(ϵ::T1, μ::T2) where {T1<:Number, T2<:Number}
    return sqrt(μ/ϵ)
end

# Gets quick indices for returning particular quadrants.
function getQuadrantSlices(numHarmonics::Integer)
    _1 = UnitRange(1,2*numHarmonics)
    _2 = UnitRange(2*numHarmonics+1,4*numHarmonics)
    return _1, _2
end
function getHalfsizeQuadrantSlices(numHarmonics::Integer)
    _1 = UnitRange(1,numHarmonics)
    _2 = UnitRange(numHarmonics+1,2*numHarmonics)
    return _1, _2
end

# Gets quick indices for returning particular quadrants.
# For a given array
function getQuadrantSlices(arr::Array{T,2}) where T<:Number
    halfsize = half( size(arr)[1] )
    _1 = UnitRange(1,halfsize)
    _2 = UnitRange(halfsize+1,2*halfsize)
    return _1, _2
end


function printFull(arr::Array{T,2}) where T<:Any
    for iy = 1:size(arr,1)
        for ix = 1:size(arr,2)
            print(arr[ix,iy], '\t')
        end
        print("\r\n")
    end
end
function printFull(filename::String, arr::Array{T,2}) where T<:Complex
    file = open(filename,"w")
    for iy = 1:size(arr,1)
        for ix = 1:size(arr,2)
            Printf.@printf(file, "'%.2f + %.2fi,", real(arr[ix,iy]), imag(arr[ix,iy]) )
            print(file,'\t')
        end
        print(file,"\r\n")
    end
end

function isEigendecomposition(A::Array{T1,2}, eigenvalues::Array{T2,2}, λ:: Array{T3,2}) where {T1<:Number, T2<:Number, T3<:Number}
    return
end


function getkz(kXY::_2VectorFloat, n::Number, wavenumber::Wavenumber, kzPositive::Bool )::ComplexF64
    nk₀ = getk₀(wavenumber)*n
    kz² = nk₀^2 - kXY[X]^2 - kXY[Y]^2
    return sqrt( kz² )*bool2posNeg(kzPositive)
end
getkz(kXY::Vector{T}, n::Number, wavenumber::Wavenumber, kzPositive::Bool ) where T<:Real = getkz(_2VectorFloat(kXY), n, wavenumber, kzPositive )

function kXYtokXYZ(kXY::_2VectorFloat, n::Number, wavenumber::Wavenumber, kzPositive::Bool )::_3VectorComplex
    kz = getkz(kXY, n, wavenumber, kzPositive )
    k₀ = getk₀(wavenumber)
    return _3VectorComplex(kXY[X], kXY[Y], kz)
end
kXYtokXYZ(kXY::Vector{T}, n::Number, wavenumber::Wavenumber, kzPositive::Bool ) where T<:Real = kXYtokXYZ(_2VectorFloat(kXY), n, wavenumber, kzPositive )


# Adds whitespace to the string so that it is the set length.  Clips the string if it is less than target length.
function setLength(string::String, len::Integer)
    if length(string) > len
        return string[1:len]
    elseif length(string) < len
        return string * " "^(len-length(string))
    end
    return string
end

function unifyStringLengths(strings::Vector{String})
    maxLength = maximum( [length(string) for string in strings] )
    unifiedStrings = [ setLength(string, maxLength) for string in strings]
    return unifiedStrings
end

function prettyPrint(arr::Array{T,2}) where T<:Number
    println("Size = ", size(arr))
    for iLine = 1:size(arr,1)
        println(arr[iLine,:])
    end
end

#Returns the fractional position of x1 between the minimum and maximum values of the input vector.  If above 1, return 1.  If below 0, return 0.
function fractionOfRange(x1::Real, minMaxX::Vector{<:Real})::Float64
    minX = minimum(minMaxX)
    maxX = maximum(minMaxX)

    if x1 > maxX
        return 1.0
    elseif x1 < minX
        return 0.0
    end
    return (x1-minX) / (maxX-minX)
end
