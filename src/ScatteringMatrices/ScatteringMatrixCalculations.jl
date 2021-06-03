# Scattering matrix calculations

# TODO: Calculate the W₀ and V₀ only once!

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
    # λᵢ::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)
    λᵢ::Diagonal{Complex{PrecisionType}} # (nHarmonics*2, nHarmonics*2)

    Vᵢ::MagneticEigenvectors{PrecisionType} # (nHarmonics*2, nHarmonics*2)

    A::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)
    B::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)
    X::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)

    S::Array{Complex{PrecisionType},2} # (nHarmonics*4, nHarmonics*4)

    _1::UnitRange{Int64}
    _2::UnitRange{Int64}
    _1Big::UnitRange{Int64}
    _2Big::UnitRange{Int64}

    function ScatteringMatrixAllocations{PrecisionType}(nHarmonics::Integer, kVectorSet::KVectorSet) where {PrecisionType<:Real}

        W₀ = ElectricEigenvectors{PrecisionType}(Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) ))
        V₀ = MagneticEigenvectors{PrecisionType}(Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) ))

        Cϵᵢⱼ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        Cμᵢⱼ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        Cϵᵢⱼ⁻¹ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        Cμᵢⱼ⁻¹ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        P = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        # P = Array{Complex{PrecisionType},2}(undef, (nHarmonics, nHarmonics) )
        Q = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        Ω² = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        Wᵢ = ElectricEigenvectors{PrecisionType}(zeros(Complex{PrecisionType}, (2*nHarmonics, 2*nHarmonics)))
        # λ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        λ = Diagonal{Complex{PrecisionType}}( 1:2*nHarmonics )

        Vᵢ = MagneticEigenvectors{PrecisionType}(zeros(Complex{PrecisionType}, (2*nHarmonics, 2*nHarmonics)))

        A = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        B = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        X = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        S = Array{Complex{PrecisionType},2}(undef, (4*nHarmonics, 4*nHarmonics) )

        _1, _2 = getQuadrantSlices(P)
        _1Big, _2Big = getQuadrantSlices(S)

        # return new(W₀, V₀, Cϵᵢⱼ, Cμᵢⱼ, Cϵᵢⱼ⁻¹, Cμᵢⱼ⁻¹, P, Q, Ω², W, λ, V, A, B, X, S)
        return new(W₀, V₀, Cϵᵢⱼ, Cμᵢⱼ, Cϵᵢⱼ⁻¹, Cμᵢⱼ⁻¹, P, Q, Ω², Wᵢ, λ, Vᵢ, A, B, X, S, _1, _2, _1Big, _2Big)
    end
end




# Calculation of P-matrix for patterned layers
# Pᵢ₁₁ = Kx_mnᵢ * Cϵ⁻¹ * Ky_mnᵢ
# Pᵢ₁₂ = Cμ - Kx_mnᵢ * Cϵ⁻¹ * Kx_mnᵢ
# Pᵢ₂₁ = Ky_mnᵢ * Cϵ⁻¹ * Ky_mnᵢ - Cμ
# Pᵢ₂₂ = -Ky_mnᵢ * Cϵ⁻¹ * Kx_mnᵢ
# Pᵢ = cat( cat(Pᵢ₁₁, Pᵢ₁₂, dims=2), cat(Pᵢ₂₁, Pᵢ₂₂, dims=2), dims=1)

# Returns the P matrix used in the eigenvalue problem for the given K-matrix and convolution matrix
function calcPmatrixPatterned( prealloc::ScatteringMatrixAllocations{PrecisionType}, kVectorSet::KVectorSet, Cϵ::Array{ComplexF64,2}, Cϵ⁻¹::Array{ComplexF64,2}, Cμ::Array{ComplexF64,2}, Cμ⁻¹::Array{ComplexF64,2}  ) where {PrecisionType<:Real}

    nHarmonics = numHarmonics(kVectorSet)

    Kx, Ky = getNormalizedKxKy(kVectorSet)

    # Shorthand for different matrix quadrants
    _1, _2 = getHalfsizeQuadrantSlices(nHarmonics)

    Pᵢ = prealloc.P
    Pᵢ[_1, _1] = Kx * Cϵ⁻¹ * Ky
    Pᵢ[_1, _2] = Cμ - Kx * Cϵ⁻¹ * Kx
    Pᵢ[_2, _1] = Ky * Cϵ⁻¹ * Ky - Cμ
    Pᵢ[_2, _2] =-Ky * Cϵ⁻¹ * Kx

    return Pᵢ
end

# Calculation of Q-matrix for patterned layers
# Qᵢ₁₁ = Kx_mnᵢ * Cμ⁻¹ * Ky_mnᵢ
# Qᵢ₁₂ = -Kx_mnᵢ * Cμ⁻¹ * Kx_mnᵢ + Cϵ
# Qᵢ₂₁ = Ky_mnᵢ * Cμ⁻¹ * Ky_mnᵢ - Cϵ
# Qᵢ₂₂ = -Ky_mnᵢ * Cμ⁻¹ * Kx_mnᵢ
# Qᵢ = cat( cat(Qᵢ₁₁, Qᵢ₁₂, dims=2), cat(Qᵢ₂₁, Qᵢ₂₂, dims=2), dims=1)
# Returns the Q matrix used in the eigenvalue problem for the given K-matrix and convolution matrix
function calcQmatrixPatterned( prealloc::ScatteringMatrixAllocations{PrecisionType}, kVectorSet::KVectorSet, Cϵ::Array{ComplexF64,2}, Cϵ⁻¹::Array{ComplexF64,2}, Cμ::Array{ComplexF64,2}, Cμ⁻¹::Array{ComplexF64,2}  ) where {PrecisionType<:Real}

    nHarmonics = numHarmonics(kVectorSet)
    Kx, Ky = getNormalizedKxKy(kVectorSet)

    # Shorthand for different matrix quadrants
    _1, _2 = getHalfsizeQuadrantSlices(nHarmonics)

    Qᵢ = prealloc.Q
    Qᵢ[_1, _1] = Kx * Cμ⁻¹ * Ky
    Qᵢ[_1, _2] =-Kx * Cμ⁻¹ * Kx + Cϵ
    Qᵢ[_2, _1] = Ky * Cμ⁻¹ * Ky - Cϵ
    Qᵢ[_2, _2] =-Ky * Cμ⁻¹ * Kx
    return Qᵢ
end


function calcΩ²( prealloc::ScatteringMatrixAllocations{PrecisionType}, P::Array{<:Number,2}, Q::Array{<:Number,2}) where {PrecisionType<:Real}

    mul!(prealloc.Ω², P, Q)  # P*Q, placing the result in prealloc.Ω²
    # Ω² = P * Q
    return prealloc.Ω²
end


# Compute eigenmodes of gap medium
# Kz = ( sqrt(I-Kx^2-Kz^2))*
# Q = [ KxKy  I-Kx^2;
#     Ky^2-I  -KxKy]
# W₀ = [I 0;
#     0 I]
# λ = [ jKz   0
#     0     jKz]
# V₀ = Q*\lambda^-1
function calcW₀(prealloc::ScatteringMatrixAllocations{PrecisionType}, numHarmonics::Integer) where {PrecisionType<:Real}
    return ElectricEigenvectors{PrecisionType}( Array{Float64,2}(I,2*numHarmonics, 2*numHarmonics) )
end
function calcW₀(numHarmonics::Integer)
    return ElectricEigenvectors{Float64}( Array{Float64,2}(I,2*numHarmonics, 2*numHarmonics) )
end

function calcWᵢλᵢ(prealloc::ScatteringMatrixAllocations{PrecisionType}, Ω²ᵢ::Array{<:Number, 2}) where {PrecisionType<:Real}
    # error()
    λ²ᵢ, Wᵢ = eigen(Ω²ᵢ)
    # @show typeof()
    # λ²ᵢ, Wᵢ = eigen!(Ω²ᵢ)  # Ω²ᵢ is no longer needed, so overwrite
    # λ²ᵢ = Diagonal( λ²ᵢ )
    # λᵢ = Diagonal( sqrt.(λ²ᵢ) )
    prealloc.λᵢ = Diagonal( sqrt.(λ²ᵢ) )
    # @show typeof(prealloc.λᵢ)
    # λᵢ = sqrt.(λ²ᵢ)
    # return ElectricEigenvectors{PrecisionType}(Wᵢ), Array(λᵢ)
    prealloc.Wᵢ = ElectricEigenvectors{PrecisionType}(Wᵢ)
    # return ElectricEigenvectors{PrecisionType}(Wᵢ), λᵢ
    return prealloc.Wᵢ, prealloc.λᵢ
end


# function Base.inv(A::Diagonal)
#     return Diagonal( inv.(A.diag) )
# end

# TODO RENAME
# function calcMagneticEigenvectorsFromQWλ(prealloc::ScatteringMatrixAllocations{PrecisionType}, Q::Array{T1,2}, Wᵢ::ElectricEigenvectors, λ::Array{T3,2}) where {T1<:Number, T2<:Number, T3<:Number, PrecisionType<:Real}
function calcMagneticEigenvectorsFromQWλ(prealloc::ScatteringMatrixAllocations{PrecisionType}, Q::AbstractArray{<:Complex,2}, Wᵢ::ElectricEigenvectors, λ::AbstractArray{<:Complex,2}) where {PrecisionType<:Real}
    # return MagneticEigenvectors{PrecisionType}( Q * Wᵢ.matrix * inv(λ) )
    # return MagneticEigenvectors{PrecisionType}( Q * Wᵢ.matrix / λ )
    # prealloc.Vᵢ.matrix = Q * Wᵢ.matrix / λ
    # @show typeof(prealloc.Vᵢ.matrix)

    # mul!(prealloc.Vᵢ.matrix, Q, Wᵢ.matrix) # Q*Wᵢ
    # @show typeof(prealloc.Vᵢ.matrix)
    # @show typeof(λ)
    # @show inv(λ)
    # rdiv!(prealloc.Vᵢ.matrix, λ) # (Q*Wᵢ) / λ
    # rdiv!(prealloc.Vᵢ.matrix, factorize(λ)) # (Q*Wᵢ) / λ
    # rdiv!(factorize(prealloc.Vᵢ.matrix), factorize(λ)) # (Q*Wᵢ) / λ
    mul!(prealloc.Vᵢ.matrix, Q, Wᵢ.matrix*inv(λ) ) # Q*(Wᵢ/λ)
    # @show typeof(prealloc.Vᵢ.matrix)
    # return MagneticEigenvectors{PrecisionType}( Q * Wᵢ.matrix / λ )
    return prealloc.Vᵢ
end

# TODO: IS THIS USED?
function calcEigenmodesForUniformLayer(prealloc::ScatteringMatrixAllocations{PrecisionType}, kVectorSet::KVectorSet, layerDef::Tlayer, matCol::MaterialCollection) where { Tlayer<:Union{UniformLayerDefinition, SemiInfiniteLayerDefinition}, PrecisionType<:Real}

    P, Q = calcPQmatrix(prealloc, layerDef, kVectorSet, matCol)
    Ω² = calcΩ²(prealloc, P, Q)
    W, λ = calcWᵢλᵢ(prealloc, Ω²)
    V₀ = calcMagneticEigenvectorsFromQWλ(prealloc, Q, W, λ)

    return V₀
end


function calcPmatrixUnpatterned(prealloc::ScatteringMatrixAllocations{PrecisionType}, kVectorSet::KVectorSet, ϵ::Number, μ::Number) where {PrecisionType<:Real}
    Kx, Ky = getNormalizedKxKy(kVectorSet)
    # @show size(Kx)
    # @show size(prealloc.P)
    # P₁₁ = Kx*Ky
    # P₁₂ = μ*ϵ*I - Kx^2
    # P₂₁ = Ky^2 - μ*ϵ*I
    # P₂₂ = -Kx*Ky
    prealloc.P[prealloc._1,prealloc._1] = (Kx*Ky) ./ ϵ
    prealloc.P[prealloc._1,prealloc._2] = (μ*ϵ*I - Kx^2) ./ ϵ
    prealloc.P[prealloc._2,prealloc._1] = (Ky^2 - μ*ϵ*I) ./ ϵ
    prealloc.P[prealloc._2,prealloc._2] = (-Kx*Ky) ./ ϵ
    # return [P₁₁ P₁₂;
    #         P₂₁ P₂₂] * (1/ϵ)
    return prealloc.P
end

# function calcQmatrixUnpatterned(prealloc::ScatteringMatrixAllocations{PrecisionType}, P::Array{ComplexF64,2}, ϵ::ComplexF64, μ::ComplexF64) where {PrecisionType<:Real}
function calcQmatrixUnpatterned(prealloc::ScatteringMatrixAllocations{PrecisionType}, P::Array{<:Complex,2}, ϵ::Complex, μ::Complex) where {PrecisionType<:Real}
    # Q = (ϵ / μ) * P
    prealloc.Q = (ϵ / μ) * P
    # return Q
    return prealloc.Q
end


# Calculate P, Q for uniform layers
function calcPQmatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::Tlayer, kVectorSet::KVectorSet,  matCol::MaterialCollection) where { Tlayer<:Union{UniformLayerDefinition, SemiInfiniteLayerDefinition}, PrecisionType<:Real}
    # ϵ, μ = calc_ϵμ(layer, matCol, kVectorSet)
    ϵ, μ = calc_ϵμ(layer, matCol, kVectorSet.wavenumber)

    ϵ = conj.(ϵ) # Not sure why the conjugate is necessary
    μ = conj.(μ)
    # lecture 7B:
    # ϵ = ϵ

    # P = calcPmatrixUnpatterned(prealloc, kVectorSet, ϵ, μ )
    P = calcPmatrixUnpatterned(prealloc, kVectorSet, ϵ, μ )
    Q = calcQmatrixUnpatterned(prealloc, P,ϵ,μ)
    return P, Q
end

function calcAB(prealloc::ScatteringMatrixAllocations{PrecisionType}, Wᵢ::ElectricEigenvectors, W₀::ElectricEigenvectors, Vᵢ::MagneticEigenvectors, V₀::MagneticEigenvectors) where {PrecisionType<:Real}
    Wᵢ⁻¹W₀ = inv(Wᵢ.matrix)*W₀.matrix
    Vᵢ⁻¹V₀ = inv(Vᵢ.matrix)*V₀.matrix
    return Wᵢ⁻¹W₀ + Vᵢ⁻¹V₀,
        Wᵢ⁻¹W₀ - Vᵢ⁻¹V₀
end


# function calcX(prealloc::ScatteringMatrixAllocations{PrecisionType}, λᵢ::Array{<:Number,2}, wavenumber::Wavenumber, thickness::Real) where {PrecisionType<:Real}
function calcX(prealloc::ScatteringMatrixAllocations{PrecisionType}, λᵢ::AbstractArray{<:Number,2}, wavenumber::Wavenumber, thickness::Real) where {PrecisionType<:Real}
    k₀d = -1*getk₀(wavenumber) * thickness
    λᵢk₀d = λᵢ*k₀d
    return exp(-λᵢ*getk₀(wavenumber)*thickness)
end


# A = W₀⁻¹Wref + V₀⁻¹Vref
function calcA_SemiInfinite(prealloc::ScatteringMatrixAllocations{PrecisionType}, W::ElectricEigenvectors, W₀::ElectricEigenvectors, V::MagneticEigenvectors, V₀::MagneticEigenvectors) where {PrecisionType<:Real}
    return inv(W₀.matrix)*W.matrix + inv(V₀.matrix)*V.matrix
end

# B = W₀⁻¹Wref - V₀⁻¹Vref
function calcB_SemiInfinite(prealloc::ScatteringMatrixAllocations{PrecisionType}, W::ElectricEigenvectors, W₀::ElectricEigenvectors, V::MagneticEigenvectors, V₀::MagneticEigenvectors) where {PrecisionType<:Real}
    return inv(W₀.matrix)*W.matrix - inv(V₀.matrix)*V.matrix
end


function calcAB_SemiInfinite(prealloc::ScatteringMatrixAllocations{PrecisionType}, W::ElectricEigenvectors, W₀::ElectricEigenvectors, V::MagneticEigenvectors, V₀::MagneticEigenvectors) where {PrecisionType<:Real}
    return calcA_SemiInfinite(prealloc, W, W₀, V, V₀), calcB_SemiInfinite(prealloc, W, W₀, V, V₀)
end


# Calculates the layer scattering matrix based on the given A, B, X values
# Aᵢ₀ = Wᵢ⁻¹W₀ + Vᵢ⁻¹V₀
# Bᵢ₀ = Wᵢ⁻¹W₀ - Vᵢ⁻¹V₀
# Xᵢ = exp(-λᵢk₀L₁)
# S₁₁ = inv(Aᵢ₀-Xᵢ*Bᵢ₀*inv(Aᵢ₀)*Xᵢ*Bᵢ₀) * (Xᵢ*Bᵢ₀*inv(Aᵢ₀)*Xᵢ*Aᵢ₀-Bᵢ₀)
# S₁₂ = inv(Aᵢ₀-Xᵢ*Bᵢ₀*inv(Aᵢ₀)*Xᵢ*Bᵢ₀) * Xᵢ*(Aᵢ₀-Bᵢ₀*inv(Aᵢ₀)*Bᵢ₀)
# S₂₁ = S₁₂
# S₂₂ = S₁₁
# function calcScatteringMatrix_ABX(prealloc::ScatteringMatrixAllocations{PrecisionType}, Aᵢ::Array{<:Number,2}, Bᵢ::Array{<:Number,2}, Xᵢ::Array{<:Number,2}) where {PrecisionType<:Real}
function calcScatteringMatrix_ABX(prealloc::ScatteringMatrixAllocations{PrecisionType}, Aᵢ::AbstractArray{<:Number,2}, Bᵢ::AbstractArray{<:Number,2}, Xᵢ::AbstractArray{<:Number,2}) where {PrecisionType<:Real}

    nHarmonics = half( size(Aᵢ)[1] )
    _1, _2 = getQuadrantSlices(nHarmonics)

    S = Array{Complex{PrecisionType},2}(undef,(4*nHarmonics,4*nHarmonics))

    if all( [ a == 0 for a in Aᵢ ] )
    # If any row is all zeros,
        @warn "calcScatteringMatrix_ABX.  One row is all zeros"
        zeroArr = zeros(Complex{PrecisionType}, (2*nHarmonics,2*nHarmonics))
        S[_1,_1] = zeroArr
        S[_1,_2] = conj(Xᵢ)
        S[_2,_1] = conj(Xᵢ)
        S[_2,_2] = zeroArr
    else
        BA⁻¹ = Bᵢ/Aᵢ
        XBA⁻¹X = Xᵢ*BA⁻¹*Xᵢ
        invAmXBA⁻¹XB = inv(Aᵢ - XBA⁻¹X*Bᵢ)

        S[_1,_1] = invAmXBA⁻¹XB * (XBA⁻¹X*Aᵢ - Bᵢ)
        S[_1,_2] = invAmXBA⁻¹XB * Xᵢ * (Aᵢ - BA⁻¹*Bᵢ)
        S[_2,_1] = S[_1,_2]
        S[_2,_2] = S[_1,_1]

    end

    return LayerScatteringMatrix(S)
end


# Calculates scattering matrix for a uniform layer
function calcScatteringMatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::UniformLayerDefinition, matCol::MaterialCollection, kVectorSet::KVectorSet) where {PrecisionType<:Real}

    W₀ = calcW₀( numHarmonics(kVectorSet) )
    V₀ = calcV₀( kVectorSet )

    P, Q = calcPQmatrix(prealloc, layer, kVectorSet, matCol)
    Ω² = calcΩ²(prealloc, P, Q)

    W, λ = calcWᵢλᵢ(prealloc,  Ω²) # Note that
    V = calcMagneticEigenvectorsFromQWλ(prealloc, Q, W, λ)

    A, B = calcAB(prealloc, W, W₀, V, V₀)
    X = calcX(prealloc, λ, kVectorSet.wavenumber, layer.thickness)

    S = calcScatteringMatrix_ABX(prealloc, A, B, X)
    return S
end


# Calculates scattering matrix for a patterned layer
function calcScatteringMatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::PatternedLayerDefinition, simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters ) where {PrecisionType<:Real}
    kVectorSet = derivedParameters.kVectorSet

    W₀ = calcW₀( numHarmonics(kVectorSet) )
    V₀ = calcV₀( kVectorSet )

    # Cϵᵢⱼ, Cμᵢⱼ = calcConvolutionMatrices( layer, lattice, gVectorSet, harmonicsSet, matCol, kVectorSet.wavenumber )
    Cϵᵢⱼ, Cμᵢⱼ = calcConvolutionMatrices( layer, simulationDefinition, derivedParameters )
    Cϵᵢⱼ⁻¹ = inv(Cϵᵢⱼ)
    Cμᵢⱼ⁻¹ = inv(Cμᵢⱼ)

    P = calcPmatrixPatterned(prealloc, kVectorSet, Cϵᵢⱼ, Cϵᵢⱼ⁻¹, Cμᵢⱼ, Cμᵢⱼ⁻¹)
    Q = calcQmatrixPatterned(prealloc, kVectorSet, Cϵᵢⱼ, Cϵᵢⱼ⁻¹, Cμᵢⱼ, Cμᵢⱼ⁻¹)
    Ω² = calcΩ²(prealloc, P, Q)

    W, λ = calcWᵢλᵢ(prealloc, Ω²)
    V = calcMagneticEigenvectorsFromQWλ(prealloc, Q, W, λ)

    A, B = calcAB(prealloc, W, W₀, V, V₀)
    X = calcX(prealloc, λ, kVectorSet.wavenumber, layer.thickness)

    S = calcScatteringMatrix_ABX(prealloc, A, B, X)
    return S
end


# Calculates the reflection layer scattering matrix based on the given A, B values
# Aᵢ₀ = W₀⁻¹Wᵢ + V₀⁻¹Vᵢ
# Bᵢ₀ = W₀⁻¹Wᵢ - V₀⁻¹Vᵢ
# S₁₁ = -Aᵢ₀⁻¹Bᵢ₀
# S₁₂ = 2*Aᵢ₀⁻¹
# S₂₁ = 0.5*(Aᵢ₀ - Bᵢ₀Aᵢ₀⁻¹Bᵢ₀)
# S₂₂ = Bᵢ₀Aᵢ₀⁻¹
function calcScatteringMatrixBottom_AB(prealloc::ScatteringMatrixAllocations{PrecisionType1}, Aᵢ₀::Array{Complex{PrecisionType2},2}, Bᵢ₀::Array{Complex{PrecisionType3},2}) where {PrecisionType1<:Real, PrecisionType2<:Real, PrecisionType3<:Real}
    PrecisionTypeNew = promote_type(PrecisionType1, PrecisionType2, PrecisionType3)

    nHarmonics = half( size(Aᵢ₀)[1] )
    _1, _2 = getQuadrantSlices(nHarmonics)

    Aᵢ₀⁻¹ = inv(Aᵢ₀)
    Bᵢ₀Aᵢ₀⁻¹ =  Bᵢ₀*Aᵢ₀⁻¹

    S = Array{Complex{PrecisionTypeNew},2}(undef,(4*nHarmonics,4*nHarmonics))

    S[_1,_1] = -Aᵢ₀⁻¹*Bᵢ₀
    S[_1,_2] = 2*Aᵢ₀⁻¹
    S[_2,_1] = 0.5*(Aᵢ₀ - Bᵢ₀Aᵢ₀⁻¹*Bᵢ₀)
    S[_2,_2] = Bᵢ₀Aᵢ₀⁻¹

    return LayerScatteringMatrix{PrecisionTypeNew}(S)
end

# Calculates the transmission layer scattering matrix based on the given A, B values
# Aᵢ₀ = Wᵢ⁻¹W₀ + Vᵢ⁻¹V₀
# Bᵢ₀ = Wᵢ⁻¹W₀ - Vᵢ⁻¹V₀
# S₁₁ = Bᵢ₀Aᵢ₀⁻¹
# S₁₂ = 0.5*(Aᵢ₀ - Bᵢ₀Aᵢ₀⁻¹Bᵢ₀)
# S₂₁ = 2*Aᵢ₀⁻¹
# S₂₂ = -Aᵢ₀⁻¹Bᵢ₀
# function calcScatteringMatrixTop_AB(prealloc::ScatteringMatrixAllocations{PrecisionType}, Aᵢ₀::Array{ComplexF64,2}, Bᵢ₀::Array{ComplexF64,2}) where {PrecisionType<:Real}
function calcScatteringMatrixTop_AB(prealloc::ScatteringMatrixAllocations{PrecisionType}, Aᵢ₀::Array{<:Complex,2}, Bᵢ₀::Array{<:Complex,2}) where {PrecisionType<:Real}

    nHarmonics = half( size(Aᵢ₀)[1] )
    _1, _2 = getQuadrantSlices(nHarmonics)

    Aᵢ₀⁻¹ = inv(Aᵢ₀)
    Bᵢ₀Aᵢ₀⁻¹ =  Bᵢ₀*Aᵢ₀⁻¹

    S = Array{ComplexF64,2}(undef,(4*nHarmonics,4*nHarmonics))

    S[_1,_1] = Bᵢ₀Aᵢ₀⁻¹
    S[_1,_2] = 0.5*(Aᵢ₀ - Bᵢ₀Aᵢ₀⁻¹*Bᵢ₀)
    S[_2,_1] = 2*Aᵢ₀⁻¹
    S[_2,_2] = -Aᵢ₀⁻¹*Bᵢ₀

    return LayerScatteringMatrix(S)
end


function calcABsemiInfiniteBottom(prealloc::ScatteringMatrixAllocations{PrecisionType}, derivedParameters::DerivedParameters, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection) where {PrecisionType<:Real}
    # @show prealloc.W₀.matrix

    # W₀ = calcW₀( numHarmonics(kVectorSet) )
    # V₀ = calcV₀( kVectorSet )
    # @show W₀
    # @show prealloc.W₀
    # @test W₀.matrix == prealloc.W₀.matrix

    P, Q = calcPQmatrix(prealloc, layer, derivedParameters.kVectorSet, matCol)
    Ω² = calcΩ²(prealloc, P, Q)
    # kz = Array(Diagonal( calckzBottom(kVectorSet, layer, matCol, kVectorSet.wavenumber) ))
    kz = Diagonal( derivedParameters.kzNormBottom )

    λ = calcΛsemiInfiniteBottom(prealloc, kz, derivedParameters.kVectorSet.wavenumber)
    # @show typeof(λ)
    # V = calcMagneticEigenvectorsFromQWλ(prealloc, Q, W₀,λ)
    V = calcMagneticEigenvectorsFromQWλ(prealloc, Q, prealloc.W₀,λ)

    A, B = calcAB_SemiInfinite(prealloc, prealloc.W₀, prealloc.W₀, V, prealloc.V₀)
    # A, B = calcAB_SemiInfinite(prealloc, prealloc.W₀, prealloc.W₀, V, prealloc.V₀)
    return A, B
end
function calcABsemiInfiniteTop(prealloc::ScatteringMatrixAllocations{PrecisionType}, derivedParameters::DerivedParameters, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection) where {PrecisionType<:Real}

    # W₀ = calcW₀( numHarmonics(kVectorSet) )
    # V₀ = calcV₀( kVectorSet )

    P, Q = calcPQmatrix(prealloc, layer, derivedParameters.kVectorSet, matCol)
    Ω² = calcΩ²(prealloc, P, Q)
    # kz = Array(Diagonal( calckzTop(kVectorSet, layer, matCol, derivedParameters.kVectorSet.wavenumber) ) )
    kz = Diagonal( derivedParameters.kzNormTop )

    λ = calcΛsemiInfiniteTop(prealloc, kz, derivedParameters.kVectorSet.wavenumber)
    V = calcMagneticEigenvectorsFromQWλ(prealloc, Q, prealloc.W₀,λ)

    A, B = calcAB_SemiInfinite(prealloc, prealloc.W₀, prealloc.W₀, V, prealloc.V₀)
    return A, B
end

# Calculates scattering matrix for the reflective layer
function calcScatteringMatrixBottom(prealloc::ScatteringMatrixAllocations{PrecisionType}, derivedParameters::DerivedParameters, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection) where {PrecisionType<:Real}
    # @show prealloc.W₀.matrix
    # A, B = calcABsemiInfiniteBottom(prealloc, layer, matCol, kVectorSet)
    # A, B = calcABsemiInfiniteBottom(prealloc, layer, matCol, kVectorSet)
    A, B = calcABsemiInfiniteBottom(prealloc, derivedParameters, layer, matCol)
    S = calcScatteringMatrixBottom_AB(prealloc, A,B)
    return S
end

# Calculates scattering matrix for a transmissive layer
# function calcScatteringMatrixTop(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection, kVectorSet::KVectorSet) where {PrecisionType<:Real}
function calcScatteringMatrixTop(prealloc::ScatteringMatrixAllocations{PrecisionType}, derivedParameters::DerivedParameters, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection) where {PrecisionType<:Real}

    # A, B = calcABsemiInfiniteTop(prealloc, layer, matCol, kVectorSet)
    A, B = calcABsemiInfiniteTop(prealloc, derivedParameters, layer, matCol)
    S = calcScatteringMatrixTop_AB(prealloc, A,B)
    return S
end

# You don't need to know the GvectorSet or lattice to calculate the scattering matrix of a uniform layer.
function calcScatteringMatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::UniformLayerDefinition, simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters ) where {PrecisionType<:Real}
    return calcScatteringMatrix(prealloc, layer, simulationDefinition.materialCollection, derivedParameters.kVectorSet )
end





# Convert the field modes into a field (2D)
# Wᵦ is NxN matrix
# modeCoeff is 2xN matrix
# Can be Wᵦ or Wₜ
function modeCoeff2Field(modeCoeff::Array{T1,2}, Wᵦ::Array{T2,2})::Array{ComplexF64,2} where {T1<:Number, T2<:Number}
    return Wᵦ*modeCoeff
end


function propagateModeCoeff(Sglobal::Array{ComplexF64,2}, sourceCoeff::Array{ComplexF64,2}, source2Coeff::Array{ComplexF64,2})

    numFields = size(source1Coeff,1)

    sourcesCoeff = vcat(source1Coeff, source2Coeff)

    scatteredCoeff = Sglobal*sourcesCoeff

    scattered1Coeff = scatteredCoeff[1:numFields,:]
    scattered2Coeff = scatteredCoeff[(numFields+1):(numFields*2),:]

    return scattered1Coeff, scattered2Coeff
end


# function calcΛsemiInfiniteBottom(prealloc::ScatteringMatrixAllocations{PrecisionType}, kz::Array{T,2}, wavenumber::Wavenumber) where {T<:Number, PrecisionType<:Real}
function calcΛsemiInfiniteBottom(prealloc::ScatteringMatrixAllocations{PrecisionType}, kz::AbstractArray{<:Complex,2}, wavenumber::Wavenumber) where {PrecisionType<:Real}
    # return vcat( hcat(-1im*kz, zeros(Complex{PrecisionType},size(kz)) ),
    #            hcat(zeros(Complex{PrecisionType},size(kz)), -1im*kz) )  # Compatible with Lecture 7

    # @show typeof(kz)
    # @test Diagonal(diag(kz)) == kz
    # @show typeof( cat(-1im*diag(kz),-1im*diag(kz); dims=1) )
    prealloc.λᵢ = Diagonal( cat(-1im*diag(kz),-1im*diag(kz),dims=1)  )
    # prealloc.λᵢ[prealloc._1,prealloc._1] = -1im*kz
    # prealloc.λᵢ[prealloc._1,prealloc._2] = zeros(Complex{PrecisionType},size(kz))
    # prealloc.λᵢ[prealloc._2,prealloc._1] = zeros(Complex{PrecisionType},size(kz))
    # prealloc.λᵢ[prealloc._2,prealloc._2] = -1im*kz

    return prealloc.λᵢ

end


# function calcΛsemiInfiniteTop(prealloc::ScatteringMatrixAllocations{PrecisionType}, kz::Array{T,2}, wavenumber::Wavenumber) where {T<:Number, PrecisionType<:Real}
function calcΛsemiInfiniteTop(prealloc::ScatteringMatrixAllocations{PrecisionType}, kz::AbstractArray{<:Complex,2}, wavenumber::Wavenumber) where {PrecisionType<:Real}
    # @show kz
    # return vcat( hcat(1im*kz, zeros(Complex{PrecisionType},size(kz)) ),
    #            hcat(zeros(Complex{PrecisionType},size(kz)), 1im*kz) )  # compatible with Lecture 7
    prealloc.λᵢ = Diagonal( cat(1im*diag(kz),1im*diag(kz),dims=1)  )
    return prealloc.λᵢ
end
# calcΛsemiInfinite(kz::LinearAlgebra.Diagonal{Complex{Float64},Array{Complex{Float64},1}}, wavenumber::Wavenumber) = calcΛsemiInfinite(Array(kz), wavenumber)


function calcV₀(prealloc::ScatteringMatrixAllocations{PrecisionType}, Q::Array{ComplexF64,2}, Λ::Array{ComplexF64,2}) where {PrecisionType<:Real}
    return MagneticEigenvectors{PrecisionType}( Q * inv(Λ) )
end
function calcV₀(Q::Array{ComplexF64,2}, Λ::Array{ComplexF64,2})
    return MagneticEigenvectors{Float64}( Q * inv(Λ) )
end


function calcV₀(kVectorSet::KVectorSet) ::MagneticEigenvectors

    KzNorm = calcFreeSpaceKzNorm(kVectorSet)
    Q = calcFreeSpaceQ( kVectorSet )
    W₀ = calcW₀( numHarmonics(kVectorSet) )
    Λ = calcFreeSpaceΛ(KzNorm)
    V₀ = calcV₀(Q,Λ)

    return V₀
end

function calcV₀(prealloc::ScatteringMatrixAllocations{PrecisionType}, PkVectorSet::KVectorSet) ::MagneticEigenvectors where {PrecisionType<:Real}

    KzNorm = calcFreeSpaceKzNorm(kVectorSet)
    Q = calcFreeSpaceQ( kVectorSet )
    W₀ = calcW₀( numHarmonics(kVectorSet) )
    Λ = calcFreeSpaceΛ(KzNorm)
    V₀ = calcV₀(prealloc,Q,Λ)

    return V₀
end
