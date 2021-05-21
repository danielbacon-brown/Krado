
# Describe the kvectors (direction for each mode/harmonic of light)
# Only one of these for a given experiment

mutable struct KVectorSet
    # Vacuum wavenumber
    wavenumber::Wavenumber

    # Array of k-vectors in the plane of the substrate. x,y coordinates.  Must be real.
    kᵢNorm::Vector{_2VectorFloat}

    # Array of k-vectors in x and y dimensions in the form of a diagonal matrix
    KxNorm::Array{Float64,2}
    KyNorm::Array{Float64,2}
end

# Define the set by listing k-vectors.  Construct the diagonal matrices
function KVectorSet(wavenumber::Wavenumber, kᵢNorm::Vector{_2VectorFloat})
    KxNorm = Diagonal( [kᵢNorm[i][XDIM] for i in 1:length(kᵢNorm) ] )
    KyNorm = Diagonal( [kᵢNorm[i][YDIM] for i in 1:length(kᵢNorm) ] )
    return KVectorSet(wavenumber, kᵢNorm, KxNorm, KyNorm)
end

# Define the set according to a central k-vector corresponding to the given harmonic order and a GvectorSet
# Each k-vector of each harmonic is just input k-vector plus the G-vector corresponding to the mode difference between that harmonic and the k-vector
# The Main way to define a k-vector set
function createKVectorSet(wavenumber::Wavenumber, k, ϖ, Gvectors::GvectorSet)
    k₀ = getk₀(wavenumber)
    k = _2VectorFloat(k)
    ϖ = _2VectorInt(ϖ)

    kᵢ = Vector{_2VectorFloat}(undef, numGvectors(Gvectors) )
    for i in 1:numGvectors(Gvectors)
        kᵢ[i] = k/k₀ + Gvectors.ΔGᵢⱼ[ Gvectors.harmonicsSet.mnᵢ[i] - ϖ ]/k₀
    end

    return KVectorSet(wavenumber, kᵢ)
end


# Define according to a 3-vector for k.
function createKVectorSet(wavenumber::Wavenumber, k::_3VectorFloat, ϖ::_2VectorInt, Gvectors::GvectorSet )
    return createKVectorSet(wavenumber, k[X:Y], ϖ, Gvectors)
end

#TODO: use defaults:

function createKVectorSet(wavenumber::Wavenumber, k, Gvectors::GvectorSet )
    return createKVectorSet(wavenumber, convert(_2VectorFloat,k[X:Y]), _2VectorInt(0,0), Gvectors)
end





# Return relative components of Kx and Ky components of K arrays
function getNormalizedKxKy(kVectorSet::KVectorSet)
    return kVectorSet.KxNorm, kVectorSet.KyNorm
end

# Return number of harmonics used
function numHarmonics(kVectorSet::KVectorSet)
    return length(kVectorSet.kᵢNorm)
end

function getk₀(kVectorSet::KVectorSet)
    return getk₀(kVectorSet.wavenumber)
end

# Returns k 2-vector of the given order index
function getkXYnorm(kVectorSet::KVectorSet, orderIndex)::_2VectorFloat
    return kVectorSet.kᵢNorm[orderIndex]
end
