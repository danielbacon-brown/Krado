# Scattering matrix calculations

# TODO: Calculate the W₀ and V₀ only once!

# Contains the data for all intermediate calculations in creating a layer scattering matrix.  Allows reuse of preallocated of data in this series of calculations
mutable struct ScatteringMatrixAllocations{PrecisionType}

    # W₀::Array{Complex{PrecisionType},2}  # size=(2*numHarmonics, 2*numHarmonics)
    # V₀::Array{Complex{PrecisionType},2}  # size=(2*numHarmonics, 2*numHarmonics)
    W₀::ElectricEigenvectors{PrecisionType}  # size=(2*numHarmonics, 2*numHarmonics)
    V₀::MagneticEigenvectors{PrecisionType}  # size=(2*numHarmonics, 2*numHarmonics)

    Cϵᵢⱼ::Array{Complex{PrecisionType},2}  # Size = (nHarmonics, nHarmonics)
    Cμᵢⱼ::Array{Complex{PrecisionType},2}  # Size = (nHarmonics, nHarmonics)
    Cϵᵢⱼ⁻¹::Array{Complex{PrecisionType},2} # Size = (nHarmonics, nHarmonics)
    Cμᵢⱼ⁻¹::Array{Complex{PrecisionType},2} # Size = (nHarmonics, nHarmonics)

    P::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)
    Q::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)

    Ω²::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)

    W::ElectricEigenvectors{PrecisionType} # (nHarmonics*2, nHarmonics*2)
    λ::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)

    V::MagneticEigenvectors{PrecisionType} # (nHarmonics*2, nHarmonics*2)

    A::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)
    B::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)
    X::Array{Complex{PrecisionType},2} # (nHarmonics*2, nHarmonics*2)

    S::Array{Complex{PrecisionType},2} # (nHarmonics*4, nHarmonics*4)

    function ScatteringMatrixAllocations{PrecisionType}(nHarmonics::Integer, kVectorSet::KVectorSet) where {PrecisionType<:Real}

        # W₀ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        # V₀ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        W₀ = ElectricEigenvectors{PrecisionType}(Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) ))
        V₀ = MagneticEigenvectors{PrecisionType}(Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) ))
        # W₀ = calcW₀(PrecisionType,nHarmonics)
        # V₀ = calcV₀(PrecisionType,kVectorSet)

        Cϵᵢⱼ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        Cμᵢⱼ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        Cϵᵢⱼ⁻¹ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        Cμᵢⱼ⁻¹ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        P = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        Q = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        Ω² = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        W = ElectricEigenvectors{PrecisionType}(zeros(Complex{PrecisionType}, (2*nHarmonics, 2*nHarmonics)))
        λ = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        V = MagneticEigenvectors{PrecisionType}(zeros(Complex{PrecisionType}, (2*nHarmonics, 2*nHarmonics)))

        A = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        B = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )
        X = Array{Complex{PrecisionType},2}(undef, (2*nHarmonics, 2*nHarmonics) )

        S = Array{Complex{PrecisionType},2}(undef, (4*nHarmonics, 4*nHarmonics) )

        return new(W₀, V₀, Cϵᵢⱼ, Cμᵢⱼ, Cϵᵢⱼ⁻¹, Cμᵢⱼ⁻¹, P, Q, Ω², W, λ, V, A, B, X, S)
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

    # Pᵢ = Array{ComplexF64,2}(undef, nHarmonics*2, nHarmonics*2)
    Pᵢ = prealloc.P
    # mul!(view(P,_1,_1)
    # KxCϵ⁻¹ = Kx * Cϵ⁻¹
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

    # Qᵢ = Array{ComplexF64,2}(undef, nHarmonics*2, nHarmonics*2)
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


# THIS FUNCTION NEVER GETS CALLED
# Calculates the electric field eigenmodes for a patterned layer
# function computeEigenModes(prealloc::ScatteringMatrixAllocations{PrecisionType}, Pᵢ::Array{ComplexF64,2}, Qᵢ::Array{ComplexF64,2}) where {PrecisionType<:Real}
#     Ω²ᵢ = calcΩ²(Pᵢ, Qᵢ)
#
#     # EMpossible: [W,LAM] = eig(OMEGA2);  LAM = sqrt(LAM);
#     # Matlab: [V,D] = eig(A) : D is diagonal matrix of eigenvalues.  V is matrix of eigenvectors.  A*V=V*D
#     # Python: w,v = eig(arr) : w is eigenvalues as vector.  v is eigenmodes as vector of vectors
#     # Julia: vals,vecs = eigen(arr)
#     # Trave11er has v_vec, Wmat = eig(Ω^2).
#     # Trave11er:
#     #v_vec, W_mat = np.linalg.eig(Omega_sq_mat) # v_vec is a vector of eigenvalues
#     #lambda_mat = np.diag(np.conj(np.sqrt(v_vec)))  # lambda_mat is a diagonal array
#     #V_mat = matmul(Q_mat, W_mat, np.linalg.inv(lambda_mat))
#     # Looks like Wᵢ=W_mat and Vᵢ=V_mat, λᵢ=lambda_mat
#
#     # conj(A) = conj.(A)
#
#     # This is the Trave11er way:
#     v_vec, W_mat = eigen(Ω²ᵢ) # v is eigenvalues as a vector
#     lambda_mat = Diagonal(conj(sqrt(v_vec))) # lambda_mat is a diagonal array from eigenvalues
#     V_mat = Qᵢ * W_mat * inv(lambda_mat)
#     return W_mat, V_mat,
#
#     # This is way I've been doing it:
#     # Calculate eigenvectors (W) and values (λ²)
#     # (Wᵢ_linear, λ²ᵢ) = eigen(Ω²ᵢ)
#     # Wᵢ = Diagonal(Wᵢ_linear)
#     # λᵢ = sqrt( λ²ᵢ )  # Way I've been doing it.  The way it's done in EMpossible
#     ### λᵢ = conj(sqrt( λ²ᵢ ))  # Adding a conjugate  # DOESN'T APPEAR TO MAKE A DIFFERENCE
#     ### λᵢ = conj.(sqrt( λ²ᵢ ))  # Adding a conjugate  # DOESN'T APPEAR TO MAKE A DIFFERENCE
#     # λᵢ⁻¹ = inv( λᵢ )
#
#     Vᵢ = Qᵢ * Wᵢ * λᵢ⁻¹
#
#     return Wᵢ, Vᵢ, λᵢ, λᵢ⁻¹
# end


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
    λ²ᵢ = Diagonal( λ²ᵢ )
    λᵢ = sqrt.(conj(λ²ᵢ)) # Based on Trave11er's RCWA.  # old
    # TODO:
    # λᵢ = sqrt.(λ²ᵢ) # EMpossible doesn't say to include a conjugate.  #TODO
    return ElectricEigenvectors{PrecisionType}(Wᵢ), Array(λᵢ)
end

# TODO RENAME
function calcMagneticEigenvectorsFromQWλ(prealloc::ScatteringMatrixAllocations{PrecisionType}, Q::Array{T1,2}, Wᵢ::ElectricEigenvectors, λ::Array{T3,2}) where {T1<:Number, T2<:Number, T3<:Number, PrecisionType<:Real}
    return MagneticEigenvectors{PrecisionType}( Q * Wᵢ.matrix * inv(λ) )
end


function calcEigenmodesForUniformLayer(prealloc::ScatteringMatrixAllocations{PrecisionType}, kVectorSet::KVectorSet, layerDef::Tlayer, matCol::MaterialCollection) where { Tlayer<:Union{UniformLayerDefinition, SemiInfiniteLayerDefinition}, PrecisionType<:Real}

    P, Q = calcPQmatrix(prealloc, layerDef, kVectorSet, matCol)
    Ω² = calcΩ²(prealloc, P, Q)
    W, λ = calcWᵢλᵢ(prealloc, Ω²)
    V₀ = calcMagneticEigenvectorsFromQWλ(prealloc, Q, W, λ)

    return V₀
end


function calcPmatrixUnpatterned(prealloc::ScatteringMatrixAllocations{PrecisionType}, kVectorSet::KVectorSet, ϵ::Number, μ::Number) where {PrecisionType<:Real}
    Kx, Ky = getNormalizedKxKy(kVectorSet)
    P₁₁ = Kx*Ky
    P₁₂ = μ*ϵ*I - Kx^2
    P₂₁ = Ky^2 - μ*ϵ*I
    P₂₂ = -Kx*Ky
    return [P₁₁ P₁₂;
            P₂₁ P₂₂] * (1/ϵ)
end

function calcQmatrixUnpatterned(prealloc::ScatteringMatrixAllocations{PrecisionType}, P::Array{ComplexF64,2}, ϵ::ComplexF64, μ::ComplexF64) where {PrecisionType<:Real}
    Q = (ϵ / μ) * P
    return Q
end


# Calculate P, Q for uniform layers
function calcPQmatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::Tlayer, kVectorSet::KVectorSet,  matCol::MaterialCollection) where { Tlayer<:Union{UniformLayerDefinition, SemiInfiniteLayerDefinition}, PrecisionType<:Real}
    ϵ, μ = calc_ϵμ(layer, matCol, kVectorSet)
    P = calcPmatrixUnpatterned(prealloc, kVectorSet, ϵ, μ )
    Q = calcQmatrixUnpatterned(prealloc, P,ϵ,μ)
    return P, Q
end

# Should be deprecated, but used for tests.
function calcA(prealloc::ScatteringMatrixAllocations{PrecisionType}, Wᵢ::ElectricEigenvectors, W₀::ElectricEigenvectors, Vᵢ::MagneticEigenvectors, V₀::MagneticEigenvectors) where {PrecisionType<:Real}
    return inv(Wᵢ.matrix)*W₀.matrix + inv(Vᵢ.matrix)*V₀.matrix
    # return inv(Wᵢ.matrix)*W₀.matrix + inv(Vᵢ.matrix)*V₀.matrix
end
function calcB(prealloc::ScatteringMatrixAllocations{PrecisionType}, Wᵢ::ElectricEigenvectors, W₀::ElectricEigenvectors, Vᵢ::MagneticEigenvectors, V₀::MagneticEigenvectors) where {PrecisionType<:Real}
    return inv(Wᵢ.matrix)*W₀.matrix - inv(Vᵢ.matrix)*V₀.matrix
end

function calcAB(prealloc::ScatteringMatrixAllocations{PrecisionType}, Wᵢ::ElectricEigenvectors, W₀::ElectricEigenvectors, Vᵢ::MagneticEigenvectors, V₀::MagneticEigenvectors) where {PrecisionType<:Real}
    Wᵢ⁻¹W₀ = inv(Wᵢ.matrix)*W₀.matrix
    Vᵢ⁻¹V₀ = inv(Vᵢ.matrix)*V₀.matrix
    return Wᵢ⁻¹W₀ + Vᵢ⁻¹V₀,
        Wᵢ⁻¹W₀ - Vᵢ⁻¹V₀
end


function calcX(prealloc::ScatteringMatrixAllocations{PrecisionType}, λᵢ::Array{<:Number,2}, wavenumber::Wavenumber, thickness::Real) where {PrecisionType<:Real}
    k₀d = -1*getk₀(wavenumber) * thickness
    λᵢk₀d = λᵢ*k₀d
    return exp(λᵢk₀d)
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
function calcScatteringMatrix_ABX(prealloc::ScatteringMatrixAllocations{PrecisionType}, Aᵢ::Array{<:Number,2}, Bᵢ::Array{<:Number,2}, Xᵢ::Array{<:Number,2}) where {PrecisionType<:Real}

    nHarmonics = half( size(Aᵢ)[1] )
    _1, _2 = getQuadrantSlices(nHarmonics)

    # Bᵢ₀Aᵢ₀⁻¹ =  Aᵢ₀\Bᵢ₀
    # invAᵢ₀mXᵢBᵢ₀Aᵢ₀XᵢBᵢ₀ = inv(Aᵢ₀ - Xᵢ*Bᵢ₀Aᵢ₀⁻¹*Xᵢ*Bᵢ₀)

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
        # BA⁻¹ = Bᵢ*inv(Aᵢ)
        BA⁻¹ = Bᵢ/Aᵢ
        XBA⁻¹X = Xᵢ*BA⁻¹*Xᵢ
        invAmXBA⁻¹XB = inv(Aᵢ - XBA⁻¹X*Bᵢ)

        S[_1,_1] = invAmXBA⁻¹XB * (XBA⁻¹X*Aᵢ - Bᵢ)
        # S[_1,_1] = inv(Aᵢ - Xᵢ*Bᵢ*inv(Aᵢ)*Xᵢ*Bᵢ) * (Xᵢ*Bᵢ*inv(Aᵢ)*Xᵢ*Aᵢ - Bᵢ)
        S[_1,_2] = invAmXBA⁻¹XB * Xᵢ * (Aᵢ - BA⁻¹*Bᵢ)
        # S[_1,_2] = inv(Aᵢ - Xᵢ*Bᵢ*inv(Aᵢ)*Xᵢ*Bᵢ) * Xᵢ * (Aᵢ - Bᵢ*inv(Aᵢ)*Bᵢ)
        S[_2,_1] = S[_1,_2]
        S[_2,_2] = S[_1,_1]

    end

    return LayerScatteringMatrix(S)
end


# More efficient construction of S:
# function calcScatteringMatrix_ABX(Aᵢ₀::Array{ComplexF64,2}, Bᵢ₀::Array{ComplexF64,2}, Xᵢ::Array{ComplexF64,2})::Array{ComplexF64,2}
#
#     # nHarmonics = round(Int64,size(Aᵢ₀)[1]/2)
#     nHarmonics = half( size(Aᵢ₀)[1] )
#     _1, _2 = getQuadrantSlices(nHarmonics)
#
#     Bᵢ₀Aᵢ₀⁻¹ =  Aᵢ₀\Bᵢ₀
#     invAᵢ₀mXᵢBᵢ₀Aᵢ₀XᵢBᵢ₀ = inv(Aᵢ₀ - Xᵢ*Bᵢ₀Aᵢ₀⁻¹*Xᵢ*Bᵢ₀)
#
#     S = Array{ComplexF64,2}(undef,(4*nHarmonics,4*nHarmonics))
#
#     S[_1,_1] = invAᵢ₀mXᵢBᵢ₀Aᵢ₀XᵢBᵢ₀ * (Xᵢ*Bᵢ₀Aᵢ₀⁻¹*Xᵢ*Aᵢ₀-Bᵢ₀)
#     S[_1,_2] = invAᵢ₀mXᵢBᵢ₀Aᵢ₀XᵢBᵢ₀ * Xᵢ*(Aᵢ₀-Bᵢ₀Aᵢ₀⁻¹*Bᵢ₀)
#     S[_2,_1] = S[_1,_2]
#     S[_2,_2] = S[_1,_1]
#
#     return S
# end

# Calculates scattering matrix for a uniform layer
function calcScatteringMatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::UniformLayerDefinition, matCol::MaterialCollection, kVectorSet::KVectorSet) where {PrecisionType<:Real}

    # Get material data for layer
    mat = getMaterial(matCol, layer.backgroundMaterialName)
    ϵ, μ = calc_ϵμ(mat, kVectorSet.wavenumber)

    W₀ = calcW₀( numHarmonics(kVectorSet) )
    V₀ = calcV₀( kVectorSet )

    P, Q = calcPQmatrix(prealloc, layer, kVectorSet, matCol)
    Ω² = calcΩ²(prealloc, P, Q)

    W, λ = calcWᵢλᵢ(prealloc,  Ω²)
    V = calcEigenmodesForUniformLayer(prealloc, kVectorSet, layer, matCol)

    A, B = calcAB(prealloc, W, W₀, V, V₀)
    X = calcX(prealloc, λ, kVectorSet.wavenumber, layer.thickness)

    S = calcScatteringMatrix_ABX(prealloc, A, B, X)
    return S

end


# Calculates scattering matrix for a patterned layer
function calcScatteringMatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::PatternedLayerDefinition, matCol::MaterialCollection, kVectorSet::KVectorSet, gVectorSet::GvectorSet, lattice::Lattice ) where {PrecisionType<:Real}

    W₀ = calcW₀( numHarmonics(kVectorSet) )
    V₀ = calcV₀( kVectorSet )

    Cϵᵢⱼ, Cμᵢⱼ = calcConvolutionMatrices( layer, lattice, gVectorSet, matCol, kVectorSet.wavenumber )
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

# Calculates scattering matrix for a uniform layer
function calcScatteringMatrixTimed(layer::PatternedLayerDefinition, matCol::MaterialCollection, kVectorSet::KVectorSet, gVectorSet::GvectorSet, lattice::Lattice )

    W₀ = calcW₀( numHarmonics(kVectorSet) )
    V₀ = calcV₀( kVectorSet )

    Cϵᵢⱼ, Cμᵢⱼ = calcConvolutionMatrices( layer, lattice, gVectorSet, matCol, kVectorSet.wavenumber )
    Cϵᵢⱼ⁻¹ = inv(Cϵᵢⱼ)
    Cμᵢⱼ⁻¹ = inv(Cμᵢⱼ)

    P = calcPmatrixPatterned(kVectorSet, Cϵᵢⱼ, Cϵᵢⱼ⁻¹, Cμᵢⱼ, Cμᵢⱼ⁻¹)
    Q = calcQmatrixPatterned(kVectorSet, Cϵᵢⱼ, Cϵᵢⱼ⁻¹, Cμᵢⱼ, Cμᵢⱼ⁻¹)

    Ω² = calcΩ²(P,Q)

    W, λ = calcWᵢλᵢ(Ω²)

    V = calcMagneticEigenvectorsFromQWλ(Q,W,λ)

    println("calc A, B")
    @time A, B = calcAB(W,W₀,V,V₀)

    println("calcX")
    # @time X = calcXold(λ, kVectorSet.wavenumber, layer.thickness)
    @time X = calcX(λ, kVectorSet.wavenumber, layer.thickness)

    println("calcS")
    @time S = calcScatteringMatrix_ABX(A,B,X)


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
    ## Bᵢ₀Aᵢ₀⁻¹ =  Aᵢ₀\Bᵢ₀
    # Bᵢ₀Aᵢ₀⁻¹ =  Aᵢ₀⁻¹*Bᵢ₀
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
function calcScatteringMatrixTop_AB(prealloc::ScatteringMatrixAllocations{PrecisionType}, Aᵢ₀::Array{ComplexF64,2}, Bᵢ₀::Array{ComplexF64,2}) where {PrecisionType<:Real}

    nHarmonics = half( size(Aᵢ₀)[1] )
    _1, _2 = getQuadrantSlices(nHarmonics)

    Aᵢ₀⁻¹ = inv(Aᵢ₀)
    # Bᵢ₀Aᵢ₀⁻¹ =  Aᵢ₀\Bᵢ₀
    # Bᵢ₀Aᵢ₀⁻¹ =  Aᵢ₀⁻¹*Bᵢ₀
    Bᵢ₀Aᵢ₀⁻¹ =  Bᵢ₀*Aᵢ₀⁻¹

    S = Array{ComplexF64,2}(undef,(4*nHarmonics,4*nHarmonics))

    S[_1,_1] = Bᵢ₀Aᵢ₀⁻¹
    S[_1,_2] = 0.5*(Aᵢ₀ - Bᵢ₀Aᵢ₀⁻¹*Bᵢ₀)
    S[_2,_1] = 2*Aᵢ₀⁻¹
    S[_2,_2] = -Aᵢ₀⁻¹*Bᵢ₀

    return LayerScatteringMatrix(S)
end

function calcABsemiInfinite(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection, kVectorSet::KVectorSet) where {PrecisionType<:Real}

    # Get material data for layer
    # λ₀ = getλ₀(kVectorSet)

    W₀ = calcW₀( numHarmonics(kVectorSet) )
    V₀ = calcV₀( kVectorSet )

    P, Q = calcPQmatrix(prealloc, layer, kVectorSet, matCol)

    Ω² = calcΩ²(prealloc, P, Q)

    # kz = Diagonal( calckz(kVectorSet, layer, matCol, kVectorSet.wavenumber) )
    kz = Array(Diagonal( calckz(kVectorSet, layer, matCol, kVectorSet.wavenumber) ))
    λ = calcΛsemiInfinite(prealloc, kz, kVectorSet.wavenumber)

    V = calcMagneticEigenvectorsFromQWλ(prealloc, Q,W₀,λ)

    A, B = calcAB_SemiInfinite(prealloc, W₀, W₀, V, V₀)
    return A, B
end


# Calculates scattering matrix for the reflective layer
function calcScatteringMatrixBottom(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection, kVectorSet::KVectorSet) where {PrecisionType<:Real}

    A, B = calcABsemiInfinite(prealloc, layer, matCol, kVectorSet)
    S = calcScatteringMatrixBottom_AB(prealloc, A,B)
    return S
end

# Calculates scattering matrix for a transmissive layer
function calcScatteringMatrixTop(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection, kVectorSet::KVectorSet) where {PrecisionType<:Real}

    A, B = calcABsemiInfinite(prealloc, layer, matCol, kVectorSet)
    S = calcScatteringMatrixTop_AB(prealloc, A,B)
    return S
end

# You don't need to know the GvectorSet or lattice to calculate the scattering matrix of a uniform layer.
function calcScatteringMatrix(prealloc::ScatteringMatrixAllocations{PrecisionType}, layer::UniformLayerDefinition, matCol::MaterialCollection, kVectorSet::KVectorSet, gVectorSet::GvectorSet, lattice::Lattice ) where {PrecisionType<:Real}
    return calcScatteringMatrix(prealloc, layer, matCol, kVectorSet )
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



function calcΛsemiInfinite(prealloc::ScatteringMatrixAllocations{PrecisionType}, kz::Array{T,2}, wavenumber::Wavenumber) where {T<:Number, PrecisionType<:Real}
    return vcat( hcat(-1im*kz, zeros(Complex{PrecisionType},size(kz)) ),
               hcat(zeros(Complex{PrecisionType},size(kz)), -1im*kz) ) * -1 / getk₀(wavenumber)
    # return Array(vcat( hcat(-1im*kz, zeros(ComplexF64,size(kz)) ),
    #            hcat(zeros(ComplexF64,size(kz)), -1im*kz) ) * -1 / getk₀(wavenumber) )
end
calcΛsemiInfinite(kz::LinearAlgebra.Diagonal{Complex{Float64},Array{Complex{Float64},1}}, wavenumber::Wavenumber) = calcΛsemiInfinite(Array(kz), wavenumber)



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
