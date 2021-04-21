
import Base.*

#TODO: Test if quadarray can perform any better than a typical array.

# Returns TRUE if the array is square (same number of rows and columns)
function isSquare(A::Array{<:Any,2})::Bool
    (m,n) = size(A)
    return m==n
end


mutable struct QuadArray{T}
    a₁₁::Array{T,2}
    a₁₂::Array{T,2}
    a₂₁::Array{T,2}
    a₂₂::Array{T,2}

    function QuadArray(B₁₁::Array{T,2}, B₁₂::Array{T,2}, B₂₁::Array{T,2}, B₂₂::Array{T,2}) where T<:Number
        # Ensure the incoming arrays are square and of same size
        @assert(size(B₁₁) == size(B₁₂) == size(B₂₁) == size(B₂₂))
        @assert(isSquare(B₁₁) && isSquare(B₁₂) && isSquare(B₂₁) && isSquare(B₂₂))
        return new{T}(B₁₁, B₁₂, B₂₁, B₂₂)
    end
end

QuadArray(B₁₁::Array{T,2}, B₁₂::Array{W,2}, B₂₁::Array{X,2}, B₂₂::Array{Y,2}) where {T,W,X,Y} = QuadArray( promote(B₁₁, B₁₂, B₂₁, B₂₂)...)


# Function to find the type that a container contains
typeOfElement(ex::QuadArray{T}) where {T} = T


function convert_eltype(::Array{T,2}, A::QuadArray{T})::Array{T,2} where T<:Number
    return [A.a₁₁ A.a₁₂; A.a₂₁ A.a₂₂ ]
end
function Array(A::QuadArray{T})::Array{T,2} where T<:Number
    return [A.a₁₁ A.a₁₂; A.a₂₁ A.a₂₂ ]
end


function *(A::QuadArray{T}, B::QuadArray{T})::QuadArray{T} where T<:Number
    AB₁₁ = A.a₁₁*B.a₁₁ + A.a₁₂*B.a₂₁
    AB₁₂ = A.a₁₁*B.a₁₂ + A.a₁₂*B.a₂₂
    AB₂₁ = A.a₂₁*B.a₁₁ + A.a₂₂*B.a₂₁
    AB₂₂ = A.a₂₁*B.a₁₂ + A.a₂₂*B.a₂₂
    AB = QuadArray(AB₁₁, AB₁₂, AB₂₁, AB₂₂)
end



function RedhefferStarProduct(A::QuadArray{T}, B::QuadArray{T})::QuadArray{T} where T<:Number
    AB₁₁ = A.a₁₁+ (A.a₁₂ * (I - inv(B.a₁₁ * A.a₂₂)) *B.a₁₁*A.a₂₁)
    AB₁₂ = A.a₁₂*inv(I - (B.a₁₁*A.a₂₂))*B.a₁₂
    AB₂₁ = B.a₁₂*inv(I - (A.a₂₂*B.a₁₁))*A.a₂₁
    AB₂₂ = B.a₂₂+ (B.a₂₁ * (I - inv(A.a₂₂ * B.a₁₁)) *A.a₂₂*B.a₁₂)
    return QuadArray(AB₁₁, AB₁₂, AB₂₁, AB₂₂)
end

# Concatenate a 2x2-array of 2d-arrays to a single 2d-array
function combineArray(arr::Array{Array{T,2},2}) where T
    return [ arr[1,1] arr[1,2]; arr[2,1] arr[2,2] ]
end

# Splits a single 2d-Array into a 2x2-array of 2d-arrays
function splitArray(arr::Array{T,2}) where T
    (lenX, lenY) = size(arr)
    @assert( iseven(lenX) && iseven(lenY))
    halfX = Int16(lenX/2)
    halfY = Int16(lenY/2)
    left = 1:halfX
    right = (halfX+1):lenX
    top = 1:halfY
    bottom = (halfY+1):lenY
    return [ [arr[top,left]] [arr[top,right]]; [arr[bottom,left]] [arr[bottom,right]] ]
end


function RedhefferStarProduct(A::Array{Array{T,2},2}, B::Array{Array{T,2},2})::Array{Array{T,2},2} where T<:Number

    AB₁₁ = A[1,1] + (A[1,2] * (I - inv(B[1,1] * A[2,2])) *B[1,1]*A[2,1])
    AB₁₂ = A[1,2]*inv(I - (B[1,1]*A[2,2]))*B[1,2]
    AB₂₁ = B[1,2]*inv(I - (A[2,2]*B[1,1]))*A[2,1]
    AB₂₂ = B[2,2] + (B[2,1] * (I - inv(A[2,2] * B[1,1])) *A[2,2]*B[1,2])
    return [ [AB₁₁] [AB₁₂]; [AB₂₁] [AB₂₂] ]
end
