# TODO: USE ADJUSTABLE TYPES
mutable struct FreeSpaceParameters

    KzNorm::Array{ComplexF64,2}   # IS THIS USED?
    Q::Array{ComplexF64,2}
    W₀::ElectricEigenvectors{Float64}
    Λ::Array{ComplexF64,2}
    V₀::MagneticEigenvectors{Float64}

    function FreeSpaceParameters(KzNorm::Array{ComplexF64,2}, Q::Array{ComplexF64,2}, W₀::ElectricEigenvectors, Λ::Array{ComplexF64,2}, V₀::MagneticEigenvectors)
        return new(KzNorm, Q, W₀, Λ, V₀)
    end

end

function FreeSpaceParameters(kVectorSet::KVectorSet)

    KzNorm = calcFreeSpaceKzNorm(kVectorSet)

    Q = calcFreeSpaceQ( kVectorSet )

    W₀ = calcW₀( numHarmonics(kVectorSet) )

    Λ = calcFreeSpaceΛ(KzNorm)  # NOTE: this is normalized

    V₀ = calcV₀(Q,Λ) # NOTE: this is normalized

    return FreeSpaceParameters(KzNorm, Q, W₀, Λ, V₀)
end



# Calculate the z-component in the layer for each k-vector
function calcFreeSpaceKzNorm(kVectorSet::KVectorSet)::Array{ComplexF64,2}
    return conj(sqrt( I -kVectorSet.KxNorm.^2 - kVectorSet.KyNorm.^2 )) # Lecture 7B
end

function calcFreeSpaceQ(kVectorSet::KVectorSet)::Array{ComplexF64,2}
    KxNorm, KyNorm = kVectorSet.KxNorm, kVectorSet.KyNorm
    return [ (KxNorm.*KyNorm) (I-KxNorm.^2);
             (KyNorm.^2-I) (-KxNorm.*KyNorm)]
end

function calcFreeSpaceΛ( KzNorm::Array{ComplexF64,2}) ::Array{ComplexF64,2}
    # DOES 'j' mean +sqrt(-1) or -sqrt(-1)?  probably plus...

    return vcat( hcat( im*KzNorm, zeros(ComplexF64,size(KzNorm)) ),
             hcat( zeros(ComplexF64,size(KzNorm)), im*KzNorm) )
end
