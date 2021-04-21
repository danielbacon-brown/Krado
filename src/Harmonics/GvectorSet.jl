# Defines the G-vectors = dot product of harmonic order and inverse lattice vector
mutable struct GvectorSet
    Gᵢ::Vector{_2VectorFloat}
    ΔGᵢⱼ::Dict{_2VectorInt,_2VectorFloat}
    harmonicsSet::HarmonicsSet
    function GvectorSet(Gᵢ::Vector{_2VectorFloat}, ΔGᵢⱼ::Dict{_2VectorInt,_2VectorFloat}, harmonicsSet::HarmonicsSet)
        return new(Gᵢ, ΔGᵢⱼ, harmonicsSet)
    end
end


function GvectorSet(Gᵢ, harmonicsSet::HarmonicsSet, lattice::Lattice)
    ΔGᵢⱼ = Dict{_2VectorInt,_2VectorFloat}()

    for Δmnᵢⱼ in harmonicsSet.Δmnᵢⱼ
        push!(ΔGᵢⱼ, Δmnᵢⱼ => Δmnᵢⱼ[1]*lattice.G₁ + Δmnᵢⱼ[2]*lattice.G₂)
    end

    return GvectorSet(Gᵢ, ΔGᵢⱼ, harmonicsSet)
end


function GvectorSet(harmonicsSet::HarmonicsSet, lattice::Lattice)
    numHarmonics = length(harmonicsSet.mnᵢ)
    Gvectors = Vector{_2VectorFloat}(undef, numHarmonics)
    Gᵢ = map( harmonic -> harmonic[U]*lattice.G₁ + harmonic[V]*lattice.G₂, harmonicsSet.mnᵢ)
    return GvectorSet(Gᵢ, harmonicsSet, lattice)
end

function numGvectors(Gvectors::GvectorSet)
    return length(Gvectors.Gᵢ)
end
function numΔGvectors(Gvectors::GvectorSet)
    return length(Gvectors.ΔGᵢⱼ)
end
