
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
