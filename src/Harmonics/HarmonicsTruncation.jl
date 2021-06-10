export HarmonicsTruncation, SuperellipseHarmonicsTruncation, γSUPERELLIPSE, calcHarmonics, HarmonicsTruncationByRectangle
# Defines how the set of harmonics to use is calculated.
abstract type HarmonicsTruncation end

# Returns harmonics set for the given truncation method
function calcHarmonicsSet(harmonicsDef::HarmonicsTruncation)
    mnᵢ = calcHarmonics(harmonicsDef)
    return HarmonicsSet( mnᵢ )
end

# Truncation method that truncates to only include harmonics (m,n) where:       (m/M)^(2*γ) + (n/N)^(2*γ) <= 1
# which enables a wide variety of shapes
mutable struct SuperellipseHarmonicsTruncation <:HarmonicsTruncation
    M::Int64
    N::Int64
    γ::Float64

    function SuperellipseHarmonicsTruncation(Mᵢ, Nᵢ, γᵢ)
        return new( Int64(Mᵢ), Int64(Nᵢ), Float64(γᵢ))
    end
end


const global γSUPERELLIPSE = Dict(
              [ ("circle", 1.0),
                ("circular", 1.0),
                ("ellipse", 1.0),
                ("pincushion", 0.3),
                ("diamond", 0.5),
                ("barrel", 1.5),
                ] )

# Calculate which harmonics meet the truncation conditions
function calcHarmonics(harmonicsDefinition::SuperellipseHarmonicsTruncation)
    M = harmonicsDefinition.M
    N = harmonicsDefinition.N
    γ = harmonicsDefinition.γ

    # Special case of only one order.
    if M == 0 && N == 0
        return [ _2VectorInt(0,0) ]
    end

    mnᵢ = _2VectorInt[]
    for n = -N:N  # inverted order
        for m = -M:M
            # (m/M)^(2*γ) + (n/N)^(2*γ) <= 1
            if ((float(m)/M)^2)^γ + ((float(n)/N)^2)^γ <= 1
                push!(mnᵢ, _2VectorInt(m,n))
            end
        end
    end
    return mnᵢ
end

# Truncation method that includes all harmonics (m,n) where -M<m<M and -N<n<N
mutable struct HarmonicsTruncationByRectangle <:HarmonicsTruncation
    M::Int64
    N::Int64
end

function calcHarmonics(harmonicsDefinition::HarmonicsTruncationByRectangle)
    M = harmonicsDefinition.M
    N = harmonicsDefinition.N

    mnᵢ = _2VectorInt[]

    for n = N:-1:-N  #inverted method
        for m = M:-1:-M
                push!(mnᵢ, _2VectorInt(m,n))
        end
    end
    return mnᵢ
end
