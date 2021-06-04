# Scattering matrix calculations






# Calculation of P-matrix for patterned layers
# Pᵢ₁₁ = Kx_mnᵢ * Cϵ⁻¹ * Ky_mnᵢ
# Pᵢ₁₂ = Cμ - Kx_mnᵢ * Cϵ⁻¹ * Kx_mnᵢ
# Pᵢ₂₁ = Ky_mnᵢ * Cϵ⁻¹ * Ky_mnᵢ - Cμ
# Pᵢ₂₂ = -Ky_mnᵢ * Cϵ⁻¹ * Kx_mnᵢ
# Pᵢ = cat( cat(Pᵢ₁₁, Pᵢ₁₂, dims=2), cat(Pᵢ₂₁, Pᵢ₂₂, dims=2), dims=1)

# Returns the P matrix used in the eigenvalue problem for the given K-matrix and convolution matrix
function calcPmatrixPatterned( prealloc::ScatteringMatrixAllocations{PrecisionType}, kVectorSet::KVectorSet, Cϵ::Array{<:Number,2}, Cϵ⁻¹::Array{<:Number,2}, Cμ::Array{<:Number,2}, Cμ⁻¹::Array{<:Number,2}  ) where {PrecisionType<:Real}
# function calcPmatrixPatterned( preallocP::AbstractArray{PrecisionType,2}, kVectorSet::KVectorSet, Cϵ::AbstractArray{<:Number,2}, Cϵ⁻¹::AbstractArray{<:Number,2}, Cμ::AbstractArray{<:Number,2}, Cμ⁻¹::AbstractArray{<:Number,2}  ) where {PrecisionType<:Real}

    # nHarmonics = numHarmonics(kVectorSet)

    Kx, Ky = getNormalizedKxKy(kVectorSet)

    prealloc.P[prealloc._1, prealloc._1] = Kx * Cϵ⁻¹ * Ky
    prealloc.P[prealloc._1, prealloc._2] = Cμ - Kx * Cϵ⁻¹ * Kx
    prealloc.P[prealloc._2, prealloc._1] = Ky * Cϵ⁻¹ * Ky - Cμ
    prealloc.P[prealloc._2, prealloc._2] =-Ky * Cϵ⁻¹ * Kx

    return prealloc.P
end

# Calculation of Q-matrix for patterned layers
# Qᵢ₁₁ = Kx_mnᵢ * Cμ⁻¹ * Ky_mnᵢ
# Qᵢ₁₂ = -Kx_mnᵢ * Cμ⁻¹ * Kx_mnᵢ + Cϵ
# Qᵢ₂₁ = Ky_mnᵢ * Cμ⁻¹ * Ky_mnᵢ - Cϵ
# Qᵢ₂₂ = -Ky_mnᵢ * Cμ⁻¹ * Kx_mnᵢ
# Qᵢ = cat( cat(Qᵢ₁₁, Qᵢ₁₂, dims=2), cat(Qᵢ₂₁, Qᵢ₂₂, dims=2), dims=1)
# Returns the Q matrix used in the eigenvalue problem for the given K-matrix and convolution matrix
function calcQmatrixPatterned( prealloc::ScatteringMatrixAllocations{PrecisionType}, kVectorSet::KVectorSet, Cϵ::AbstractArray{<:Number,2}, Cϵ⁻¹::AbstractArray{<:Number,2}, Cμ::AbstractArray{<:Number,2}, Cμ⁻¹::AbstractArray{<:Number,2}  ) where {PrecisionType<:Real}

    Kx, Ky = getNormalizedKxKy(kVectorSet)

    prealloc.Q[prealloc._1, prealloc._1] = Kx * Cμ⁻¹ * Ky
    prealloc.Q[prealloc._1, prealloc._2] =-Kx * Cμ⁻¹ * Kx + Cϵ
    prealloc.Q[prealloc._2, prealloc._1] = Ky * Cμ⁻¹ * Ky - Cϵ
    prealloc.Q[prealloc._2, prealloc._2] =-Ky * Cμ⁻¹ * Kx
    return prealloc.Q
end


function calcΩ²( prealloc::ScatteringMatrixAllocations{PrecisionType}, P::AbstractArray{<:Number,2}, Q::AbstractArray{<:Number,2}) where {PrecisionType<:Real}

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
    λ²ᵢ, Wᵢ = eigen(Ω²ᵢ)
    prealloc.λᵢ = Diagonal( sqrt.(λ²ᵢ) )
    prealloc.Wᵢ.matrix = Wᵢ
    return prealloc.Wᵢ, prealloc.λᵢ
end

# TODO RENAME
function calcMagneticEigenvectorsFromQWλ(prealloc::ScatteringMatrixAllocations{PrecisionType}, Q::AbstractArray{<:Complex,2}, Wᵢ::ElectricEigenvectors, λ::AbstractArray{<:Complex,2}) where {PrecisionType<:Real}
    mul!(prealloc.Vᵢ.matrix, Q, Wᵢ.matrix*inv(λ) ) # Q*(Wᵢ/λ)

    return prealloc.Vᵢ
end


function calcPmatrixUnpatterned(prealloc::ScatteringMatrixAllocations{PrecisionType}, kVectorSet::KVectorSet, ϵ::Number, μ::Number) where {PrecisionType<:Real}
    Kx, Ky = getNormalizedKxKy(kVectorSet)

    prealloc.P[prealloc._1,prealloc._1] = (Kx*Ky) ./ ϵ
    prealloc.P[prealloc._1,prealloc._2] = (μ*ϵ*I - Kx^2) ./ ϵ
    prealloc.P[prealloc._2,prealloc._1] = (Ky^2 - μ*ϵ*I) ./ ϵ
    prealloc.P[prealloc._2,prealloc._2] = (-Kx*Ky) ./ ϵ
    return prealloc.P
end

function calcQmatrixUnpatterned(prealloc::ScatteringMatrixAllocations{PrecisionType}, P::Array{<:Complex,2}, ϵ::Complex, μ::Complex) where {PrecisionType<:Real}
    # Q = (ϵ / μ) * P
    prealloc.Q = (ϵ / μ) * P
    return prealloc.Q
end


# Calculate P, Q for uniform layers
function calcPQmatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::Tlayer, kVectorSet::KVectorSet,  matCol::MaterialCollection) where { Tlayer<:Union{UniformLayerDefinition, SemiInfiniteLayerDefinition}, PrecisionType<:Real}
    ϵ, μ = calc_ϵμ(layer, matCol, kVectorSet.wavenumber)

    ϵ = conj.(ϵ) # Not sure why the conjugate is necessary
    μ = conj.(μ)
    # lecture 7B:
    # ϵ = ϵ

    P = calcPmatrixUnpatterned(prealloc, kVectorSet, ϵ, μ )
    Q = calcQmatrixUnpatterned(prealloc, P,ϵ,μ)
    return P, Q
end

# TODO: Make this preallocated?
function calcAB(prealloc::ScatteringMatrixAllocations{PrecisionType}, Wᵢ::ElectricEigenvectors, W₀::ElectricEigenvectors, Vᵢ::MagneticEigenvectors, V₀::MagneticEigenvectors) where {PrecisionType<:Real}
    Wᵢ⁻¹W₀ = Wᵢ.matrix \ W₀.matrix
    Vᵢ⁻¹V₀ = Vᵢ.matrix\V₀.matrix
    return Wᵢ⁻¹W₀ + Vᵢ⁻¹V₀,
           Wᵢ⁻¹W₀ - Vᵢ⁻¹V₀
end


function calcX(prealloc::ScatteringMatrixAllocations{PrecisionType}, λᵢ::AbstractArray{<:Number,2}, wavenumber::Wavenumber, thickness::Real) where {PrecisionType<:Real}
    k₀d = -1*getk₀(wavenumber) * thickness
    rmul!(λᵢ, k₀d)
    prealloc.X = exp(λᵢ)
    return prealloc.X
end


# A = W₀⁻¹Wref + V₀⁻¹Vref
# B = W₀⁻¹Wref - V₀⁻¹Vref
# For Semi-infinite: W == W₀
function calcABfromWV_SemiInfinite(prealloc::ScatteringMatrixAllocations{PrecisionType}, V::MagneticEigenvectors, V₀::MagneticEigenvectors) where {PrecisionType<:Real}
    V₀⁻¹V = inv(V₀.matrix)*V.matrix
    A = I + V₀⁻¹V
    B = I - V₀⁻¹V
    return A, B
end


# Calculates the layer scattering matrix based on the given A, B, X values
# Aᵢ₀ = Wᵢ⁻¹W₀ + Vᵢ⁻¹V₀
# Bᵢ₀ = Wᵢ⁻¹W₀ - Vᵢ⁻¹V₀
# Xᵢ = exp(-λᵢk₀L₁)
# S₁₁ = inv(Aᵢ₀-Xᵢ*Bᵢ₀*inv(Aᵢ₀)*Xᵢ*Bᵢ₀) * (Xᵢ*Bᵢ₀*inv(Aᵢ₀)*Xᵢ*Aᵢ₀-Bᵢ₀)
# S₁₂ = inv(Aᵢ₀-Xᵢ*Bᵢ₀*inv(Aᵢ₀)*Xᵢ*Bᵢ₀) * Xᵢ*(Aᵢ₀-Bᵢ₀*inv(Aᵢ₀)*Bᵢ₀)
# S₂₁ = S₁₂
# S₂₂ = S₁₁
function calcScatteringMatrix_ABX(prealloc::ScatteringMatrixAllocations{PrecisionType}, Aᵢ::AbstractArray{<:Number,2}, Bᵢ::AbstractArray{<:Number,2}, Xᵢ::AbstractArray{<:Number,2}) where {PrecisionType<:Real}

    nHarmonics = half( size(Aᵢ)[1] )
    # _1, _2 = getQuadrantSlices(nHarmonics)


    if all( [ a == 0 for a in Aᵢ ] )
    # If any row is all zeros,
        @warn "calcScatteringMatrix_ABX.  One row is all zeros"
        zeroArr = zeros(Complex{PrecisionType}, (2*nHarmonics,2*nHarmonics))
        prealloc.S.matrix[prealloc._1Big,prealloc._1Big] = zeroArr
        prealloc.S.matrix[prealloc._1Big,prealloc._2Big] = conj(Xᵢ)
        prealloc.S.matrix[prealloc._2Big,prealloc._1Big] = conj(Xᵢ)
        prealloc.S.matrix[prealloc._2Big,prealloc._2Big] = zeroArr
    else
        BA⁻¹ = Bᵢ/Aᵢ
        XBA⁻¹X = Xᵢ*BA⁻¹*Xᵢ
        invAmXBA⁻¹XB = inv(Aᵢ - XBA⁻¹X*Bᵢ)

        prealloc.S.matrix[prealloc._1Big,prealloc._1Big] = invAmXBA⁻¹XB * (XBA⁻¹X*Aᵢ - Bᵢ)
        prealloc.S.matrix[prealloc._1Big,prealloc._2Big] = invAmXBA⁻¹XB * Xᵢ * (Aᵢ - BA⁻¹*Bᵢ)
        prealloc.S.matrix[prealloc._2Big,prealloc._1Big] = prealloc.S.matrix[prealloc._1Big,prealloc._2Big]
        prealloc.S.matrix[prealloc._2Big,prealloc._2Big] = prealloc.S.matrix[prealloc._1Big,prealloc._1Big]

    end

    return prealloc.S
end


# Calculates scattering matrix for a uniform layer
function calcScatteringMatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::UniformLayerDefinition, matCol::MaterialCollection, kVectorSet::KVectorSet) where {PrecisionType<:Real}

    P, Q = calcPQmatrix(prealloc, layer, kVectorSet, matCol)
    Ω² = calcΩ²(prealloc, P, Q)

    W, λ = calcWᵢλᵢ(prealloc,  Ω²) # Note that
    V = calcMagneticEigenvectorsFromQWλ(prealloc, Q, W, λ)

    A, B = calcAB(prealloc, W, prealloc.W₀, V, prealloc.V₀)
    X = calcX(prealloc, λ, kVectorSet.wavenumber, layer.thickness)

    S = calcScatteringMatrix_ABX(prealloc, A, B, X)
    return S
end


# Calculates scattering matrix for a patterned layer
function calcScatteringMatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::PatternedLayerDefinition, simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters ) where {PrecisionType<:Real}
    kVectorSet = derivedParameters.kVectorSet

    Cϵᵢⱼ, Cμᵢⱼ = calcConvolutionMatrices( prealloc.Cϵᵢⱼ, prealloc.Cμᵢⱼ, layer, simulationDefinition, derivedParameters )
    prealloc.Cϵᵢⱼ⁻¹[:,:] = inv(Cϵᵢⱼ)
    prealloc.Cμᵢⱼ⁻¹[:,:] = inv(Cμᵢⱼ)

    P = calcPmatrixPatterned(prealloc, kVectorSet, prealloc.Cϵᵢⱼ, prealloc.Cϵᵢⱼ⁻¹, prealloc.Cμᵢⱼ, prealloc.Cμᵢⱼ⁻¹)
    Q = calcQmatrixPatterned(prealloc, kVectorSet, prealloc.Cϵᵢⱼ, prealloc.Cϵᵢⱼ⁻¹, prealloc.Cμᵢⱼ, prealloc.Cμᵢⱼ⁻¹)
    Ω² = calcΩ²(prealloc, P, Q)

    W, λ = calcWᵢλᵢ(prealloc, Ω²)
    V = calcMagneticEigenvectorsFromQWλ(prealloc, Q, W, λ)

    A, B = calcAB(prealloc, W, prealloc.W₀, V, prealloc.V₀)
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
function calcScatteringMatrixBottom_AB(prealloc::ScatteringMatrixAllocations{PrecisionType}, Aᵢ₀::AbstractArray{<:Complex,2}, Bᵢ₀::Array{<:Complex,2}) where {PrecisionType<:Real}

    nHarmonics = half( size(Aᵢ₀)[1] )

    Aᵢ₀⁻¹ = inv(Aᵢ₀)
    Bᵢ₀Aᵢ₀⁻¹ =  Bᵢ₀*Aᵢ₀⁻¹

    prealloc.S.matrix[prealloc._1Big,prealloc._1Big] = -Aᵢ₀⁻¹*Bᵢ₀
    prealloc.S.matrix[prealloc._1Big,prealloc._2Big] = 2*Aᵢ₀⁻¹
    prealloc.S.matrix[prealloc._2Big,prealloc._1Big] = 0.5*(Aᵢ₀ - Bᵢ₀Aᵢ₀⁻¹*Bᵢ₀)
    prealloc.S.matrix[prealloc._2Big,prealloc._2Big] = Bᵢ₀Aᵢ₀⁻¹

    return prealloc.S
end

# Calculates the transmission layer scattering matrix based on the given A, B values
# Aᵢ₀ = Wᵢ⁻¹W₀ + Vᵢ⁻¹V₀
# Bᵢ₀ = Wᵢ⁻¹W₀ - Vᵢ⁻¹V₀
# S₁₁ = Bᵢ₀Aᵢ₀⁻¹
# S₁₂ = 0.5*(Aᵢ₀ - Bᵢ₀Aᵢ₀⁻¹Bᵢ₀)
# S₂₁ = 2*Aᵢ₀⁻¹
# S₂₂ = -Aᵢ₀⁻¹Bᵢ₀
function calcScatteringMatrixTop_AB(prealloc::ScatteringMatrixAllocations{PrecisionType}, Aᵢ₀::Array{<:Complex,2}, Bᵢ₀::Array{<:Complex,2}) where {PrecisionType<:Real}

    Aᵢ₀⁻¹ = inv(Aᵢ₀)
    Bᵢ₀Aᵢ₀⁻¹ =  Bᵢ₀*Aᵢ₀⁻¹

    prealloc.S.matrix[prealloc._1Big,prealloc._1Big] = Bᵢ₀Aᵢ₀⁻¹
    prealloc.S.matrix[prealloc._1Big,prealloc._2Big] = 0.5*(Aᵢ₀ - Bᵢ₀Aᵢ₀⁻¹*Bᵢ₀)
    prealloc.S.matrix[prealloc._2Big,prealloc._1Big] = 2*Aᵢ₀⁻¹
    prealloc.S.matrix[prealloc._2Big,prealloc._2Big] = -Aᵢ₀⁻¹*Bᵢ₀

    return prealloc.S
end


function calcABsemiInfiniteBottom(prealloc::ScatteringMatrixAllocations{PrecisionType}, derivedParameters::DerivedParameters, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection) where {PrecisionType<:Real}

    P, Q = calcPQmatrix(prealloc, layer, derivedParameters.kVectorSet, matCol)
    Ω² = calcΩ²(prealloc, P, Q)
    kz = Diagonal( derivedParameters.kzNormBottom )

    λ = calcΛsemiInfiniteBottom(prealloc, derivedParameters.kzNormBottom)

    V = calcMagneticEigenvectorsFromQWλ(prealloc, Q, prealloc.W₀,λ)

    A, B = calcABfromWV_SemiInfinite(prealloc, V, prealloc.V₀)
    return A, B
end
function calcABsemiInfiniteTop(prealloc::ScatteringMatrixAllocations{PrecisionType}, derivedParameters::DerivedParameters, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection) where {PrecisionType<:Real}

    P, Q = calcPQmatrix(prealloc, layer, derivedParameters.kVectorSet, matCol)
    Ω² = calcΩ²(prealloc, P, Q)

    λ = calcΛsemiInfiniteTop(prealloc, derivedParameters.kzNormTop)
    V = calcMagneticEigenvectorsFromQWλ(prealloc, Q, prealloc.W₀,λ)

    A, B = calcABfromWV_SemiInfinite(prealloc, V, prealloc.V₀)
    return A, B
end

# Calculates scattering matrix for the reflective layer
function calcScatteringMatrixBottom(prealloc::ScatteringMatrixAllocations{PrecisionType}, derivedParameters::DerivedParameters, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection) where {PrecisionType<:Real}
    A, B = calcABsemiInfiniteBottom(prealloc, derivedParameters, layer, matCol)
    S = calcScatteringMatrixBottom_AB(prealloc, A,B)
    return S
end

# Calculates scattering matrix for a transmissive layer
function calcScatteringMatrixTop(prealloc::ScatteringMatrixAllocations{PrecisionType}, derivedParameters::DerivedParameters, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection) where {PrecisionType<:Real}

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


function calcΛsemiInfiniteBottom(prealloc::ScatteringMatrixAllocations{PrecisionType}, kz::AbstractArray{<:Complex,1}) where {PrecisionType<:Real}
    prealloc.λᵢ = Diagonal( cat(-1im*kz,-1im*kz,dims=1)  )
    return prealloc.λᵢ
end


function calcΛsemiInfiniteTop(prealloc::ScatteringMatrixAllocations{PrecisionType}, kz::AbstractArray{<:Complex,1}) where {PrecisionType<:Real}
    prealloc.λᵢ = Diagonal( cat(1im*kz,1im*kz,dims=1)  )
    return prealloc.λᵢ
end


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
