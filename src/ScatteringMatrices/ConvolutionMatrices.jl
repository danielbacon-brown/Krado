

# Calculates the convolution matrices ⟦ϵ⟧  ⟦μ⟧ a the given layer and lattice with the corresponding harmonics and materials
function calcConvolutionMatrices( preallocCϵᵢⱼ, preallocCμᵢⱼ, layerDef::PatternedLayerDefinition, simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters )

    positionGrid, ϵGrid, μGrid = getPositionϵμGrids(layerDef, simulationDefinition.lattice, simulationDefinition.materialCollection, getWavenumber(simulationDefinition))

    # UNSURE - For some reason taking the conjugate is necessary.  Does not affect IntegrationTest3 but is needed for absorptive films
    ϵGrid = conj.(ϵGrid)
    μGrid = conj.(μGrid)

    Cϵᵢⱼ = calcConvolutionMatrix( preallocCϵᵢⱼ, ϵGrid, positionGrid, derivedParameters)
    Cμᵢⱼ = calcConvolutionMatrix( preallocCμᵢⱼ, μGrid, positionGrid, derivedParameters)

    return Cϵᵢⱼ, Cμᵢⱼ
end

# Input: dict describing ∫ϵ⋅exp(Δ𝐆⋅𝐫)
#    and the set describing the harmonics to use
# Output: 2D grid describing the convolution matrix (e.g. ⟦ϵ⟧  or ⟦μ⟧ )
# Iterate over every combination of G-vectors and grab appropriate value from dict
function assembleConvolutionMatrix( preallocCxᵢⱼ::AbstractArray{<:Number,2}, valuesByΔϖ::Dict{_2VectorInt,ComplexF64}, harmonicsSet::HarmonicsSet  )
    # Refers to either Cϵ, ⟦ϵ⟧  or Cμ, ⟦μ⟧

    for i_G = 1:numHarmonics(harmonicsSet)
        for j_G = 1:numHarmonics(harmonicsSet)
            preallocCxᵢⱼ[i_G,j_G] = valuesByΔϖ[ harmonicsSet.mnᵢ[i_G] - harmonicsSet.mnᵢ[j_G]]
        end
    end
    return preallocCxᵢⱼ
end

# Calculates the convolution matrix for x (either ϵ or μ) by creating a dictionary of results for all possible pairs of G-vectors, then quickly generating the matrix using the dict
function calcConvolutionMatrix( preallocCxᵢⱼ, xGrid, positionGrid::PositionGridXY, derivedParameters::DerivedParameters)
    ∫xexpΔ𝐆𝐫Dict = calc∫xexpΔ𝐆𝐫Dict(xGrid, positionGrid, derivedParameters)
    preallocCxᵢⱼ = assembleConvolutionMatrix( preallocCxᵢⱼ, ∫xexpΔ𝐆𝐫Dict, derivedParameters.harmonicsSet )
    return preallocCxᵢⱼ
end


# Calculate the transform ∫∫x⋅exp(G̅⋅𝐫) for all G-vectorDifferences, putting results in a dict with the harmonic mn as the key.
# This corresponds to the 'a' matrix
# valueGrid is a 2D grid of the values at r real-space coordinates.  Does transform using a single G-vector: G.
function calc∫xexpΔ𝐆𝐫Dict(valuesGrid, positionGrid::PositionGridXY, Gvectors::GvectorSet, harmonicsSet::HarmonicsSet)

    ∫xexpΔ𝐆𝐫 = Dict{_2VectorInt, ComplexF64}()
    for ΔharmonicIndex in 1:numΔGvectors(Gvectors)
        harmonic = harmonicsSet.Δmnᵢⱼ[ΔharmonicIndex]
        ΔG = Gvectors.ΔGᵢⱼ[harmonic]
        push!( ∫xexpΔ𝐆𝐫, harmonic => calc∫xexp𝐆𝐫(valuesGrid, positionGrid, ΔG) )
    end
    return ∫xexpΔ𝐆𝐫
end
function calc∫xexpΔ𝐆𝐫Dict(valuesGrid, positionGrid::PositionGridXY, derivedParameters::DerivedParameters)
    return calc∫xexpΔ𝐆𝐫Dict(valuesGrid, positionGrid, derivedParameters.gVectorSet, derivedParameters.harmonicsSet)
end
