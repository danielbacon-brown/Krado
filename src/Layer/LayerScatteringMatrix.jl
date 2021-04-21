# Represents the scattering matrix of a particular layer.
mutable struct LayerScatteringMatrix{PrecisionType<:Real}

    # Size = (4*nHarmonics) x (4*nHarmonics)
    # matrix::A
    matrix::Array{Complex{PrecisionType},2}
    
    function LayerScatteringMatrix{PrecisionType}(matrix::Array{Complex{PrecisionType}}) where {PrecisionType<:Real}
        @assert isSquare(matrix)
        return new(matrix)
    end
    
end

function LayerScatteringMatrix(matrix::Array{Complex{PrecisionType}}) where {PrecisionType<:Real}
    @assert isSquare(matrix)
    return LayerScatteringMatrix{PrecisionType}(matrix)
end

# import Base.size
function Base.size(sm::LayerScatteringMatrix)
    return size(sm.matrix)
end

function numHarmonics(layerScatteringMatrix::LayerScatteringMatrix)
    return quarter(size(layerScatteringMatrix,1))
end

function getQuadrantSlices(sm::LayerScatteringMatrix) 
    return getQuadrantSlices(sm.matrix)
end

function Base.copy(sm::LayerScatteringMatrix)
    return LayerScatteringMatrix(copy(sm.matrix))
end

# Returns four views, each referring to a quadrant of the array
function getQuadViews(A)
    _1, _2 = getQuadrantSlices(A)
    A₁₁ = view(A, _1, _1)
    A₁₂ = view(A, _1, _2)
    A₂₁ = view(A, _2, _1)
    A₂₂ = view(A, _2, _2)
    return A₁₁, A₁₂, A₂₁, A₂₂
end



# Redheffer star product of A and B, with result stored in A.
# All inputs must be square and same size
#
# AB₁₁ = A₁₁ + (A₁₂ * inv(I - (B₁₁ * A₂₂)) *B₁₁*A₂₁)
# AB₁₂ = A₁₂*inv(I - (B₁₁*A₁₂))*B₁₂
# AB₂₁ = B₂₁*inv(I - (A[_2,_2]*B₁₁))*A₂₁
# AB₂₂ = B₂₂ + (B₂₁ * inv(I - (A₂₂ * B₁₁)) *A₂₂*B₁₂)
#
# Note that each section's calculation does not depend on a section that was previously calculated.  This means each section can be stored as an in-place operation A->AB.
function RedhefferStarProduct!( AscatteringMatrix::LayerScatteringMatrix, BscatteringMatrix::LayerScatteringMatrix)
    
    @inbounds begin    
        A₁₁, A₁₂, A₂₁, A₂₂ = getQuadViews(AscatteringMatrix.matrix)
        B₁₁, B₁₂, B₂₁, B₂₂ = getQuadViews(BscatteringMatrix.matrix)
        
        A₁₂_ImB₁₁A₂₂⁻¹ = A₁₂ * inv(I - (B₁₁ * A₂₂))  # TODO: preallocate
        B₂₁_ImA₂₂B₁₁⁻¹ = B₂₁ * inv(I - (A₂₂ * B₁₁))
        
        A₁₁[:,:] = A₁₁ .+ A₁₂_ImB₁₁A₂₂⁻¹ * B₁₁ * A₂₁
        A₁₂[:,:] = A₁₂_ImB₁₁A₂₂⁻¹ * B₁₂
        A₂₁[:,:] = B₂₁_ImA₂₂B₁₁⁻¹ * A₂₁
        A₂₂[:,:] = B₂₂ .+ B₂₁_ImA₂₂B₁₁⁻¹ * A₂₂ * B₁₂
    end;
        
    return nothing
end



# Redheffer star product of A and B, with result stored in a new array.  This is more memory intensive.
# All inputs must be square and same size
#
# AB₁₁ = A₁₁ + (A₁₂ * inv(I - (B₁₁ * A₂₂)) *B₁₁*A₂₁)
# AB₁₂ = A₁₂*inv(I - (B₁₁*A₁₂))*B₁₂
# AB₂₁ = B₂₁*inv(I - (A[_2,_2]*B₁₁))*A₂₁
# AB₂₂ = B₂₂ + (B₂₁ * inv(I - (A₂₂ * B₁₁)) *A₂₂*B₁₂)
function RedhefferStarProduct(AscatteringMatrix::LayerScatteringMatrix, BscatteringMatrix::LayerScatteringMatrix)::LayerScatteringMatrix 
    A = AscatteringMatrix.matrix
    B = BscatteringMatrix.matrix

    _1, _2 = getQuadrantSlices(A)
    AB = Array{ComplexF64,2}(undef,size(A))

    AB[_1,_1] = A[_1,_1] + (A[_1,_2] * inv(I - (B[_1,_1] * A[_2,_2])) *B[_1,_1]*A[_2,_1])
    AB[_1,_2] = A[_1,_2]*inv(I - (B[_1,_1]*A[_2,_2]))*B[_1,_2]
    AB[_2,_1] = B[_2,_1]*inv(I - (A[_2,_2]*B[_1,_1]))*A[_2,_1]
    AB[_2,_2] = B[_2,_2] + (B[_2,_1] * inv(I - (A[_2,_2] * B[_1,_1])) *A[_2,_2]*B[_1,_2])
    return LayerScatteringMatrix(AB)
end

# Syntactic sugar for
function ⊗(A::LayerScatteringMatrix, B::LayerScatteringMatrix)::LayerScatteringMatrix
    return RedhefferStarProduct(A, B)
end
function ⊗(A, B) 
    return RedhefferStarProduct(A, B)
end
