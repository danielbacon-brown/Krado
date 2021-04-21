
# Describe the kvectors (direction for each mode/harmonic of light)
# Only one of these for a given experiment
mutable struct KVectorSet
    # Vacuum wavelength and wavenumber of the incident vector
    # λ₀::Float64
    # k₀::Float64
    wavenumber::Wavenumber

    # Array of k-vectors in the plane of the substrate. x,y coordinates.  Must be real.
    kᵢ::Vector{_2VectorFloat}

    # Array of k-vectors in x and y dimensions in the form of a diagonal matrix
    Kx::Array{Float64,2}
    Ky::Array{Float64,2}
end

# Define the set by listing k-vectors.  Construct the diagonal matrices
function KVectorSet(wavenumber::Wavenumber, kᵢ::Vector{_2VectorFloat})

    Kx = Diagonal( [kᵢ[i][XDIM] for i in 1:length(kᵢ) ] )
    Ky = Diagonal( [kᵢ[i][YDIM] for i in 1:length(kᵢ) ] )

    return KVectorSet(wavenumber, kᵢ, Kx, Ky)
end

# Define the set according to a central k-vector corresponding to the given harmonic order and a GvectorSet
# Each k-vector of each harmonic is just input k-vector plus the G-vector corresponding to the mode difference between that harmonic and the k-vector
# The Main way to define a k-vector set
function createKVectorSet(wavenumber::Wavenumber, k, ϖ, Gvectors::GvectorSet)
    k = _2VectorFloat(k)
    ϖ = _2VectorInt(ϖ)

    kᵢ = Vector{_2VectorFloat}(undef, numGvectors(Gvectors) )
    for i in 1:numGvectors(Gvectors)
        kᵢ[i] = k + Gvectors.ΔGᵢⱼ[ Gvectors.harmonicsSet.mnᵢ[i] - ϖ ]
    end

    return KVectorSet(wavenumber, kᵢ)
end
# createKVectorSet(wavenumber::Wavenumber, k::TU2VectorReal, ϖ::TU2VectorInt, Gvectors::GvectorSet) = createKVectorSet(wavenumber, _2VectorFloat(k), _2VectorInt(ϖ), Gvectors)


# Define according to a 3-vector for k.
function createKVectorSet(wavenumber::Wavenumber, k::_3VectorFloat, ϖ::_2VectorInt, Gvectors::GvectorSet )
    return createKVectorSet(wavenumber, k[X:Y], ϖ, Gvectors)
end
# createKVectorSet(wavenumber::Wavenumber, k, ϖ, Gvectors::GvectorSet ) = createKVectorSet(wavenumber::Wavenumber, _3VectorFloat(k), _2VectorInt(ϖ), Gvectors::GvectorSet )

#TODO: use defaults

function createKVectorSet(wavenumber::Wavenumber, k, Gvectors::GvectorSet )
    return createKVectorSet(wavenumber, convert(_2VectorFloat,k[X:Y]), _2VectorInt(0,0), Gvectors)
end
# createKVectorSet(wavenumber::Wavenumber, k, Gvectors::GvectorSet ) = createKVectorSet(wavenumber,_3VectorFloat(k), Gvectors::GvectorSet )





# Return relative components of Kx and Ky components of K arrays
function getNormalizedKxKy(kVectorSet::KVectorSet)
    k₀ = getk₀(kVectorSet.wavenumber)
    return kVectorSet.Kx/k₀, kVectorSet.Ky/k₀
end

# Return number of harmonics used
function numHarmonics(kVectorSet::KVectorSet)
    return length(kVectorSet.kᵢ)
end

function getk₀(kVectorSet::KVectorSet)
    return getk₀(kVectorSet.wavenumber)
end

# Returns k 2-vector of the given order index
function getkXY(kVectorSet::KVectorSet, orderIndex)::_2VectorFloat
    return kVectorSet.kᵢ[orderIndex]
end
