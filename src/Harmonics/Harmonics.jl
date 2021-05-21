
# Create a set that contains (mnᵢ - mnⱼ) for all i,j in input
function createDifferenceVectors(mnᵢ)
    Δmnᵢⱼ = Set{_2VectorInt}()
    numHarmonics = length(mnᵢ)
    # Create a set, eliminating duplicates
    for i = 1:numHarmonics
        for j = 1:numHarmonics
            push!(Δmnᵢⱼ, mnᵢ[j] - mnᵢ[i])
        end
    end
    # Convert to Vector to maintain order, and sort
    Δmnᵢⱼ = sort([k for k in Δmnᵢⱼ])
    return Δmnᵢⱼ
end


# Defines a set of harmonics m,n to use in the simulation.  Relates the (m,n) harmonic to the harmonic index according to a dict
mutable struct HarmonicsSet
    # List of (m,n) tuples
    mnᵢ::Vector{_2VectorInt}
    # Dict relates the key (m,n) to the corresponding index in mnᵢ
    indᵢ_mn::Dict{_2VectorInt,Int64}

    # List of all possible difference vectors of the set mnᵢ (mnᵢ - mnⱼ)
    Δmnᵢⱼ::Vector{_2VectorInt}

    function HarmonicsSet(mnᵢ::Vector{_2VectorInt}, indᵢ_mn::Dict{_2VectorInt,Int64}, Δmnᵢⱼ::Vector{_2VectorInt})
        # Make sure the vectors and tuples are self consistent
        @assert all( [ indᵢ_mn[ mnᵢ[ind] ] == ind for ind = 1:length(mnᵢ)] )

        return new(mnᵢ, indᵢ_mn, Δmnᵢⱼ)
    end

end

function HarmonicsSet(mnᵢ::Vector{_2VectorInt}, indᵢ_mn::Dict{_2VectorInt,Int64})
    # Make sure the vectors and tuples are self consistent
    @assert all( [ indᵢ_mn[ mnᵢ[ind] ] == ind for ind = 1:length(mnᵢ)] )

    Δmnᵢ = createDifferenceVectors(mnᵢ)
    return HarmonicsSet(mnᵢ, indᵢ_mn, Δmnᵢ)
end

# Creates set using separate vectors of m, n, combining into a vector of tuples and then inverting than into a dict
function HarmonicsSet(mnᵢ)
    mnᵢ = convert(Vector{_2VectorInt},mnᵢ)
    indᵢ_mn = Dict( mnᵢ .=> 1:length(mnᵢ))
    return HarmonicsSet(mnᵢ, indᵢ_mn)
end

function numHarmonics(harmonicsSet::HarmonicsSet)
    return length(harmonicsSet.mnᵢ)
end
function numΔHarmonics(harmonicsSet::HarmonicsSet)
    return length(harmonicsSet.Δmnᵢⱼ)
end

# Returns the index of the included harmonic order
function getOrderIndex(harmonicsSet::HarmonicsSet, ϖ::_2VectorInt)
    ϖ = _2VectorInt(ϖ)

    local index
    try
        index = harmonicsSet.indᵢ_mn[ϖ]
    catch err
        ϖX, ϖY = ϖ[X], ϖ[Y]
        error("Requested index is not in harmonics set: ($ϖX,$ϖY)" )
    end
    return index
end
