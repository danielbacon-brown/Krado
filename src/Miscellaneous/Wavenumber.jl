# Represents the vacuum wavenumber/wavelength/photon energy
# Uses SI units
# f*λ₀ = c
# k₀ = 2*π/λ₀
# k₀ = E/(ħc)
# ω = 2*π*f
# T = 1/f

# Uses wavenumber as default creator to improve understandability
mutable struct Wavenumber
    k₀::Float64
    function Wavenumber(k₀::Real)
        @assert k₀ > 0
        new(Float64(k₀))
    end
end

# Creation functions
function WavenumberByk₀(k₀::Real)
    return Wavenumber(k₀)
end

function WavenumberByλ₀(λ₀::Real)
    k₀ = 2*π/λ₀
    return WavenumberByk₀(k₀)
end

function WavenumberByf(f::Real)
    λ₀ = c/f
    return WavenumberByλ₀(λ₀)
end

function WavenumberByω(ω::Real)
    f = 2*π/ω
    return WavenumberByf(f)
end

function WavenumberByT(T::Real)  # Temporal Period
    f = 1/T
    return WavenumberByf(f)
end

function WavenumberByE(E::Real)
    k₀ = E/(ħ*c)
    return WavenumberByk₀( k₀ )
end



# Conversion functions
function getk₀(wavenumber::Wavenumber)
    return wavenumber.k₀
end

function getλ₀(wavenumber::Wavenumber)
    return 2*π/getk₀(wavenumber)
end

function getf(wavenumber::Wavenumber)
    return c/getλ₀(wavenumber)
end

function getω(wavenumber::Wavenumber)
    return 2*π*getf(wavenumber)
end

function getT(wavenumber::Wavenumber)
    return 1/getf(wavenumber)
end

function getE(wavenumber::Wavenumber)
    return getk₀(wavenumber) * ħ * c
end


function Base.isapprox(x::Wavenumber, y::Wavenumber; rtol::Real=sqrt(eps), atol::Real=0, nans::Bool=false, norm::Function)
    return isapprox(x.k₀,y.k₀; rtol=rtol, atol=atol, nans=nans, norm=norm)
end
