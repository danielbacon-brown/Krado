
export BoundaryDefinition, InputByOrderBoundaryDefinition, getWavenumber

abstract type BoundaryDefinition end

mutable struct InputByOrderBoundaryDefinition <:BoundaryDefinition

    # Define the characteristics of the incident light
    # The azimuthal and zenith angles defining the direction of incidence.  The light is travelling at this angle within the substrate or superstrate material.
    wavenumber::Wavenumber
    θ::Float64 # Polar angle
    ϕ::Float64 # Azimuthal angle
    mainHarmonicOrder::_2VectorInt # Which harmonic does this represent (0,0 order by default)
    isTop::Bool # Is the light coming from the top or bottom? This only matters for calculating the central k-vector from the angles

    # Dictionary defining the mode amplitudes, A, for each integer harmonic order, ϖ
    # Amplitudes are essentially Jones vectors [Aₛ, Aₚ]
    Abyϖbottom::Dict{_2VectorInt, _2VectorComplex}
    Abyϖtop::Dict{_2VectorInt, _2VectorComplex}

    function InputByOrderBoundaryDefinition(wavenumber::Wavenumber, θ::Real, ϕ::Real, mainHarmonicOrder::_2VectorInt, isTop::Bool, Abyϖbottom::Dict{_2VectorInt,_2VectorComplex}, Abyϖtop::Dict{_2VectorInt,_2VectorComplex})

        # No longer necessary with algorithm updates:
        # If the incidence angle is exactly normal or very close to it, the A-matrix of a uniform Air layer can become singular.  Set the angle so this doesn't happen.
        # if -1e-6 < θ < 1e-6
        #     θ = (θ==0) ? 1e-6 : 1e-6*sign(θ)
        # end
        # if -1e-6 < ϕ < 1e-6
        #     ϕ = (ϕ==0) ? 1e-6 : 1e-6*sign(ϕ)
        # end

        return new(wavenumber, Float64(θ), Float64(ϕ), mainHarmonicOrder, isTop, Abyϖbottom, Abyϖtop)
    end
end
InputByOrderBoundaryDefinition(wavenumber::Wavenumber, θ::Real, ϕ::Real, mainHarmonicOrder::TU2VectorInt, isTop::Bool, Abyϖbottom::Dict{_2VectorInt,_2VectorComplex}, Abyϖtop::Dict{_2VectorInt,_2VectorComplex}) = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, convert(_2VectorInt,mainHarmonicOrder), isTop, Abyϖbottom, Abyϖtop)


function InputByOrderBoundaryDefinition(wavenumber::Wavenumber, θ::Real, ϕ::Real,  isTop::Bool, amplitudes; mainHarmonicOrder=_2VectorInt(0,0))
    if isTop
        Abyϖtop = Dict{_2VectorInt, _2VectorComplex}( mainHarmonicOrder => amplitudes )
        Abyϖbottom = Dict{_2VectorInt, _2VectorComplex}()
    else
        Abyϖbottom = Dict{_2VectorInt, _2VectorComplex}( mainHarmonicOrder => amplitudes )
        Abyϖtop = Dict{_2VectorInt, _2VectorComplex}()
    end
    return InputByOrderBoundaryDefinition(wavenumber, Float64(θ), Float64(ϕ), mainHarmonicOrder, isTop, Abyϖbottom, Abyϖtop)
end



function getWavenumber(inputByOrderBoundaryDefinition::InputByOrderBoundaryDefinition)
    return inputByOrderBoundaryDefinition.wavenumber
end
