# Define materials

# Any kind of material.  Returns (ϵ,μ) based on wavelength and material and (x,y)
abstract type AbstractMaterial end

# Return ϵ or μ based on wavelength and position
abstract type PermittivityModel end
abstract type PermeabilityModel end


#If a more particular function is not defined for this, then remove the position in call
function calc_ϵ(permittivity::PermittivityModel, wavenumber::Wavenumber, position)
    return calc_ϵ(permittivity, wavenumber)
end
# calc_ϵ(permittivity::PermittivityModel, wavenumber::Wavenumber, position) = calc_ϵ(permittivity, wavenumber, _2VectorFloat(position) )

function calc_μ(permeability::PermeabilityModel, wavenumber::Wavenumber, position)
    return calc_μ(permeability, wavenumber)
end
# calc_μ(permittivity::PermittivityModel, wavenumber::Wavenumber, position) = calc_μ(permittivity, wavenumber, _2VectorFloat(position) )



# Permittivity is allowed to be complex
struct ConstantPermittivity <: PermittivityModel
    ϵ::ComplexF64
    function ConstantPermittivity(ϵᵢ::Number)
        return new(convert(ComplexF64,ϵᵢ) )
    end
end

function calc_ϵ(permittivity::ConstantPermittivity, wavenumber::Wavenumber)
    return permittivity.ϵ
end


# Permeability is allowed to be complex
struct ConstantPermeability <: PermeabilityModel
    μ::ComplexF64
    function ConstantPermeability(μᵢ::Number)
        return new(convert(ComplexF64,μᵢ) )
    end
end
function calc_μ(permeability::ConstantPermeability, wavenumber::Wavenumber)
    return permeability.μ
end



# Defines permittivity as a function taking (wavenumber::AbstractFloat)
struct FunctionPermittivity{T<:Function} <: PermittivityModel
    ϵfunction::T
    function FunctionPermittivity{T}(ϵfunction::T) where {T<:Function}
        return new{T}(ϵfunction)
    end
end
function FunctionPermittivity(ϵfunction::T) where {T<:Function}
    return FunctionPermittivity{T}(ϵfunction)
end

function calc_ϵ(permittivity::FunctionPermittivity, wavenumber::Wavenumber)
    return permittivity.ϵfunction(wavenumber)
end


# Defines permittivity as a function taking (wavenumber::AbstractFloat, position::TU2VectorReal)
struct SpatialFunctionPermittivity{T<:Function} <: PermittivityModel
    ϵfunction::T
    function SpatialFunctionPermittivity{T}(ϵfunction::T) where {T<:Function}
        return new{T}(ϵfunction)
    end
end
function SpatialFunctionPermittivity(ϵfunction::T) where {T<:Function}
    return SpatialFunctionPermittivity{T}(ϵfunction)
end


function calc_ϵ(permittivity::SpatialFunctionPermittivity, wavenumber::Wavenumber, position)
    return permittivity.ϵfunction(wavenumber, position)
end
# calc_ϵ(permittivity::SpatialFunctionPermittivity, wavenumber::Wavenumber, position) = calc_ϵ(permittivity, wavenumber, _2VectorFloat(position))



# struct InterpolatedPermittivity <: PermittivityModel
#     interpolation::Interpolations.Extrapolation
#     function InterpolatedPermittivity(interpolation::Interpolations.Extrapolation)
#         return new(interpolation)
#     end
# end
# 
# function calc_ϵ(permittivity::InterpolatedPermittivity, wavenumber::Wavenumber)
#     itp = permittivity.interpolation
#     return itp( getλ₀(wavenumber) )
# end


const global WAVENUMBERLIMITSDEFAULT = (WavenumberByk₀(floatmax(Float64)), WavenumberByk₀(floatmin(Float64)))

mutable struct Material <: AbstractMaterial
    ϵModel::PermittivityModel
    μModel::PermeabilityModel
    wavenumberRange::Tuple{Wavenumber,Wavenumber}  # The order of the wavenumber range does not matter.
    function Material(ϵᵢ::PermittivityModel, μᵢ::PermeabilityModel = ConstantPermeability(1.0); wavenumberRange = WAVENUMBERLIMITSDEFAULT)
        @assert length(wavenumberRange) == 2
        wavenumberRange = (wavenumberRange[1], wavenumberRange[2]) # convert to tuple
        
        return new(ϵᵢ, μᵢ, wavenumberRange)
    end
end

# Error if the wavenumber is outside the given range
function checkWavenumberRange(material::Material, wavenumber::Wavenumber)
    λ₀ = getλ₀(wavenumber)
    λ₁ = getλ₀(material.wavenumberRange[1])
    λ₂ = getλ₀(material.wavenumberRange[2])
    if (λ₀ != λ₁) && (λ₀ != λ₂) && ((λ₀ < λ₁) == (λ₀ < λ₂))
        throw(DomainError(λ₀, "Wavelength $λ₀ not in valid range $λ₁ - $λ₂ m."))
    end
end

# Calls to material get sent to the permittivity/permeability model
function calc_ϵ(material::Material, wavenumber::Wavenumber)
    checkWavenumberRange(material, wavenumber)
    return calc_ϵ(material.ϵModel, wavenumber)
end
function calc_μ(material::Material, wavenumber::Wavenumber)
    checkWavenumberRange(material, wavenumber)
    return calc_μ(material.μModel, wavenumber)
end
function calc_ϵμ(material::Material, wavenumber::Wavenumber)
    return calc_ϵ(material, wavenumber), calc_μ(material, wavenumber)
end

# Collection at a single point
function calc_ϵ(material::Material, wavenumber::Wavenumber, position)
    checkWavenumberRange(material, wavenumber)
    return calc_ϵ(material.ϵModel, wavenumber, position)
end
function calc_μ(material::Material, wavenumber::Wavenumber, position)
    checkWavenumberRange(material, wavenumber)
    return calc_μ(material.μModel, wavenumber, position)
end
function calc_ϵμ(material::Material, wavenumber::Wavenumber, position)
    return calc_ϵ(material, wavenumber, position), calc_μ(material, wavenumber, position)
end
# calc_ϵ(material::Material, wavenumber::Wavenumber, position::TU2VectorReal) = calc_ϵ(material, wavenumber, _2VectorFloat(position))
# calc_μ(material::Material, wavenumber::Wavenumber, position::TU2VectorReal) = calc_μ(material, wavenumber, _2VectorFloat(position))
# calc_ϵμ(material::Material, wavenumber::Wavenumber, position::TU2VectorReal) = calc_ϵμ(material, wavenumber, _2VectorFloat(position))

# Syntactic sugar for getting n, z
# Returns the complex refractive index
function calc_n(material::AbstractMaterial, wavenumber::Wavenumber, position)
    ϵ, μ = calc_ϵμ(material, wavenumber, position)
    return convert_ϵμ2n(ϵ, μ)
end
# calc_n(material::AbstractMaterial, wavenumber::Wavenumber, position::TU2VectorReal) = calc_n(material, wavenumber, _2VectorFloat(position))

# Returns the complex refractive index
function calc_n(material::AbstractMaterial, wavenumber::Wavenumber)
    ϵ, μ = calc_ϵμ(material, wavenumber)
    return convert_ϵμ2n(ϵ, μ)
end

# Returns the complex impedance
function calc_z(material::AbstractMaterial, wavenumber::Wavenumber, position)
    ϵ, μ = calc_ϵμ(material, wavenumber, position)
    return convert_ϵμ2z(ϵ, μ)
end
# calc_z(material::AbstractMaterial, wavenumber::Wavenumber, position::TU2VectorReal) = calc_z(material, wavenumber, _2VectorFloat(position) )

# Returns the complex impedance
function calc_z(material::AbstractMaterial, wavenumber::Wavenumber)::ComplexF64
    ϵ, μ = calc_ϵμ(material, wavenumber)
    return convert_ϵμ2z(ϵ, μ)
end


# Calculates the ϵ and μ for each position in a grid of materials
# This is not very efficient, but that probably doesn't matter.
function calc_ϵμ(gridMaterials, wavenumber::Wavenumber) 
    gridϵ = map( mat -> calc_ϵ(mat, wavenumber), gridMaterials)
    gridμ = map( mat -> calc_μ(mat, wavenumber), gridMaterials)
    return gridϵ, gridμ
end


# Adjusts the k-vector by the material's refractive index.  May be complex
function k₀2kInMaterial(wavenumber::Wavenumber, mat::AbstractMaterial)
    return getk₀(wavenumber) * calc_n(mat, wavenumber )
end
