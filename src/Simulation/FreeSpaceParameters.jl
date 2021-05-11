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

    Λ = calcFreeSpaceΛ(KzNorm)

    V₀ = calcV₀(Q,Λ)

    return FreeSpaceParameters(KzNorm, Q, W₀, Λ, V₀)
end


# Calculate the z-component in the layer for each k-vector
function calcFreeSpaceKzNorm(kVectorSet::KVectorSet)::Array{ComplexF64,2}
    KxNorm, KyNorm = getNormalizedKxKy(kVectorSet)
    return conj(sqrt( I - KxNorm.^2 - KyNorm.^2 ))
    # return conj(sqrt( I - KxNorm^2 - KyNorm^2 )) # old
end

function calcFreeSpaceQ(kVectorSet::KVectorSet)::Array{ComplexF64,2}
    KxNorm, KyNorm = getNormalizedKxKy(kVectorSet)
    return [ (KxNorm*KyNorm) (I-KxNorm^2);
             (KyNorm^2-I) (-KxNorm*KyNorm)]
end

function calcFreeSpaceΛ( KzNorm::Array{ComplexF64,2}) ::Array{ComplexF64,2}
    # old
    return vcat( hcat( im*KzNorm, zeros(ComplexF64,size(KzNorm)) ),
             hcat( zeros(ComplexF64,size(KzNorm)), im*KzNorm) )

    # Lecture 7B test?
    # return vcat( hcat( 1im*conj(KzNorm), zeros(ComplexF64,size(KzNorm)) ),
    #          hcat( zeros(ComplexF64,size(KzNorm)), 1im*conj(KzNorm)) )  # same result as above.  Probably due to the complex components of kzNorm not mattering for no evanescent waves
    # return vcat( hcat( -1im*conj(KzNorm), zeros(ComplexF64,size(KzNorm)) ), #Breaks everything
    #          hcat( zeros(ComplexF64,size(KzNorm)), -1im*conj(KzNorm)) )
    # return vcat( hcat( -1im*KzNorm, zeros(ComplexF64,size(KzNorm)) ), #Breaks everything
    #          hcat( zeros(ComplexF64,size(KzNorm)), -1im*KzNorm) )
end

# function calcV₀(Q::Array{ComplexF64,2}, Λ::Array{ComplexF64,2})
#     return MagneticEigenvectors{Float64}( Q * inv(Λ) )
# end
# # MOVED
# function calcV₀(prealloc::ScatteringMatrixAllocations{PrecisionType}, Q::Array{ComplexF64,2}, Λ::Array{ComplexF64,2})
#     return MagneticEigenvectors{PrecisionType}( Q * inv(Λ) )
# end

# function calcV₀(kVectorSet::KVectorSet) ::MagneticEigenvectors
#
#     KzNorm = calcFreeSpaceKzNorm(kVectorSet)
#     Q = calcFreeSpaceQ( kVectorSet )
#     W₀ = calcW₀( numHarmonics(kVectorSet) )
#     Λ = calcFreeSpaceΛ(KzNorm)
#     V₀ = calcV₀(Q,Λ)
#
#     return V₀
# end
#
# function calcV₀(prealloc::ScatteringMatrixAllocations{PrecisionType}, PkVectorSet::KVectorSet) ::MagneticEigenvectors where {PrecisionType<:Real}
#
#     KzNorm = calcFreeSpaceKzNorm(kVectorSet)
#     Q = calcFreeSpaceQ( kVectorSet )
#     W₀ = calcW₀( numHarmonics(kVectorSet) )
#     Λ = calcFreeSpaceΛ(KzNorm)
#     V₀ = calcV₀(prealloc,Q,Λ)
#
#     return V₀
# end
