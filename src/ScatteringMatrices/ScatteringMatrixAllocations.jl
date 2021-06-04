# Contains the data for all intermediate calculations in creating a layer scattering matrix.  Allows reuse of preallocated of data in this series of calculations
mutable struct ScatteringMatrixAllocations{PrecisionType}

    W₀::ElectricEigenvectors{PrecisionType}  # size=(2*numHarmonics, 2*numHarmonics)
    V₀::MagneticEigenvectors{PrecisionType}  # size=(2*numHarmonics, 2*numHarmonics)

    Cϵᵢⱼ::Array{Complex{PrecisionType},2}  # Size = (nHarmonics, nHarmonics)
    Cμᵢⱼ::Array{Complex{PrecisionType},2}  # Size = (nHarmonics, nHarmonics)
    Cϵᵢⱼ⁻¹::Array{Complex{PrecisionType},2} # Size = (nHarmonics, nHarmonics)
    Cμᵢⱼ⁻¹::Array{Complex{PrecisionType},2} # Size = (nHarmonics, nHarmonics)

    P::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)
    Q::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)

    Ω²::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)

    Wᵢ::ElectricEigenvectors{PrecisionType} # (nHarmonics*2, nHarmonics*2)
    λᵢ::Diagonal{Complex{PrecisionType}} # (nHarmonics*2, nHarmonics*2)

    Vᵢ::MagneticEigenvectors{PrecisionType} # (nHarmonics*2, nHarmonics*2)

    A::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)
    B::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)
    X::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)

    S::LayerScatteringMatrix{PrecisionType} # (nHarmonics*4, nHarmonics*4)
    Sassembled::LayerScatteringMatrix{PrecisionType} # (nHarmonics*4, nHarmonics*4)

    _1::UnitRange{Int64}
    _2::UnitRange{Int64}
    _1Big::UnitRange{Int64}
    _2Big::UnitRange{Int64}

    function ScatteringMatrixAllocations{PrecisionType}(nHarmonics::Integer, kVectorSet::KVectorSet) where {PrecisionType<:Real}

        W₀ = ElectricEigenvectors{PrecisionType}(Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) ))
        V₀ = MagneticEigenvectors{PrecisionType}(Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) ))

        Cϵᵢⱼ = Array{Complex{PrecisionType},2}(undef, (nHarmonics, nHarmonics) )
        Cμᵢⱼ = Array{Complex{PrecisionType},2}(undef, (nHarmonics, nHarmonics) )
        Cϵᵢⱼ⁻¹ = Array{Complex{PrecisionType},2}(undef, (nHarmonics, nHarmonics) )
        Cμᵢⱼ⁻¹ = Array{Complex{PrecisionType},2}(undef, (nHarmonics, nHarmonics) )

        P = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        Q = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        Ω² = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        Wᵢ = ElectricEigenvectors{PrecisionType}(zeros(Complex{PrecisionType}, (2*nHarmonics, 2*nHarmonics)))
        λ = Diagonal{Complex{PrecisionType}}( 1:2*nHarmonics )

        Vᵢ = MagneticEigenvectors{PrecisionType}(zeros(Complex{PrecisionType}, (2*nHarmonics, 2*nHarmonics)))

        A = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        B = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        X = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        S = LayerScatteringMatrix( Array{Complex{PrecisionType},2}(undef, (4*nHarmonics, 4*nHarmonics) ) )
        Sassembled = LayerScatteringMatrix( Array{Complex{PrecisionType},2}(undef, (4*nHarmonics, 4*nHarmonics) ) )

        _1, _2 = getQuadrantSlices(P)
        _1Big, _2Big = getQuadrantSlices(S)

        return new(W₀, V₀, Cϵᵢⱼ, Cμᵢⱼ, Cϵᵢⱼ⁻¹, Cμᵢⱼ⁻¹, P, Q, Ω², Wᵢ, λ, Vᵢ, A, B, X, S, Sassembled, _1, _2, _1Big, _2Big)
    end
end
