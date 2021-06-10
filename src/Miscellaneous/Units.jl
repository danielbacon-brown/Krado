# Define Units and constants

export _2VectorInt, _2VectorFloat, _2VectorComplex, _3VectorInt, _3VectorFloat, _3VectorComplex
export TU2VectorReal
export m, μm, cm, mm, μm, nm, Å, s, A, kg, degrees, radians
export FORWARD, BACKWARD, BOTTOMINDEX, TOPINDEX, BOTTOM, TOP
export c, ϵ, ϵ₀, h, ħ


# Dimensions used in arrays
const global XDIM = 1
const global YDIM = 2
const global ZDIM = 3
const global X = 1
const global Y = 2
const global Z = 3

# Dimensions in terms of lattice vectors
const global U = 1
const global V = 2

# Dimensions in terms of ŝ, p̂
const global S = 1
const global P = 2

const global FORWARD = true
const global BACKWARD = false

const global BOTTOM = false
const global TOP = true
const global BOTTOMINDEX = 1
const global TOPINDEX = 2

const global m = 1
const global cm = 1e-2
const global mm = 1e-3
const global μm = 1e-6
const global nm = 1e-9
const global Å = 1e-10

const global LENGTHLABEL = Dict{Float64,String}(
   [ (1.0, "m"),
   (1e-2, "cm"),
   (1e-3, "mm"),
   (1e-6, "μm"),
   (1e-9, "nm"),
   (1e-10, "Å"),
   ] )

# Si base units
const global s = 1 # seconds
const global A = 1 # amps
const global kg = 1 # kilograms

# SI derived units  # SHOULD I ACTUALLY USE THESE?
const global F = (s^4 * A^2) / (m^2*kg) # Farads
const global J = kg*m^2/s^2 # Joules

const global degrees = π/180
const global radians = 1

const global MINIMUMAOILIMIT = 1e-5

# NOTE: Spherical coordinates follow the (ρ, θ, ϕ) notation described in https://mathworld.wolfram.com/SphericalCoordinates.html

const global c = 2.99792458e8 * m / s # Vacuum speed of light
const global ϵ₀ = 8.8541878176e-12 * F/m # F/m =  A2⋅s4⋅kg−1⋅m−3
const global h = 6.62607015e-34 * J*s # Planck's constant
const global ħ = 1.054571817e-34 * J*s # Planck's constant

# Plotting
# const global PATCHES = pyimport("matplotlib.patches")
@show PATCHES
@show PATCHES
const global LATTICECOLOR = [0.6, 0.6, 0.6]
const global PLOTTINGOFFSET = 0.01
const global KVECTORCOLOR = [0.0, 0.0, 0.0]
const global KVECTORREALLINESTYLE = "-"
const global KVECTORIMAGLINESTYLE = "--"
const global KVECTORWIDTH = 1e-8
const global PVECTORCOLOR = [0.0, 0.5, 0.0]
const global PVECTORREALLINESTYLE = "-"
const global PVECTORIMAGLINESTYLE = "dashed"
const global PVECTORWIDTH = 0.5e-8

# Defining syntactic sugar for small vectors types
const _2VectorInt = SVector{2,Int64}  # 2-vector containing Int64's
const _2VectorFloat = SVector{2,Float64} # 2-vector containing Float64's.
const _2VectorComplex = SVector{2,ComplexF64} # 2-vector containing Float64's
const _3VectorInt = SVector{3,Int64} # 3-vector containing ints
const _3VectorFloat = SVector{3,Float64}
const _3VectorComplex = SVector{3,ComplexF64} # 3-vector containing complex floats

# Type unions: if a function is called with nonstandard types, the call should default to a different function that converts to the standard types and yields the desired function
const TU2VectorInt = Union{Vector{T}, SVector{2,T}} where {T<:Integer}
const TU2VectorReal = Union{Vector{T}, SVector{2,T}} where {T<:Real}
const TU2VectorComplex = Union{Vector{T}, SVector{2,T}} where {T<:Number}
const TU3VectorInt = Union{Vector{T}, SVector{3,T}} where {T<:Integer}
const TU3VectorReal = Union{Vector{T}, SVector{3,T}} where {T<:Real}
const TU3VectorComplex = Union{Vector{T}, SVector{3,T}} where {T<:Number}

const TUArrayLikeComplex = Union{ Array{T1,2}, LinearAlgebra.Diagonal{T2,Array{T3,1}}} where {T1<:Number, T2<:Number, T3<:Number}

# Unit vectors
const x̂ = _3VectorInt(1,0,0)
const ŷ = _3VectorInt(0,1,0)
const ẑ = _3VectorInt(0,0,1)


const global surfaceNormal = _3VectorInt(0,0,1)


const global DEFAULTSEQUENTIALCOLORMAP = get_cmap(:plasma)

const global EPSILON = 1
const global MU = 2


const global DEFAULTPLOTSIZE = (5,5)
