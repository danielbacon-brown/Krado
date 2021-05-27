# Harmonics related plotting functions
function setLabels(ax, xLabel, yLabel)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
end

# Scatterplot showing the (m,n) harmonic orders
function plotHarmonicsSet(simulation::SimulationDefinition)

    derivedParameters = DerivedParameters(simulation)
    harmonicsSet = derivedParameters.harmonicsSet

    fig, ax = create2Dfigure(title="Harmonics (m,n)")

    nPoints = numHarmonics(harmonicsSet)
    mᵢ = Vector{Int64}(undef,nPoints)
    nᵢ = Vector{Int64}(undef,nPoints)
    for iHarmonic in 1:nPoints
        mᵢ[iHarmonic] = harmonicsSet.mnᵢ[iHarmonic][U]
        nᵢ[iHarmonic] = harmonicsSet.mnᵢ[iHarmonic][V]
    end

    ax.scatter(mᵢ, nᵢ )
    setLabels(ax, "m", "n")

    return fig, ax
end


# Scatterplot showing the values of all G-vectors
function plotGvectorSet(simulation::SimulationDefinition; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    derivedParameters = DerivedParameters(simulation)
    gVectorSet = derivedParameters.gVectorSet

    fig, ax = create2Dfigure(title="G-vectors")

    nPoints = numGvectors(gVectorSet)
    Gᵢx = Vector{Float64}(undef,nPoints)
    Gᵢy = Vector{Float64}(undef,nPoints)
    for iHarmonic in 1:nPoints
        Gᵢx[iHarmonic] = gVectorSet.Gᵢ[iHarmonic][X] * scale
        Gᵢy[iHarmonic] = gVectorSet.Gᵢ[iHarmonic][Y] * scale
    end

    ax.scatter(Gᵢx, Gᵢy )
    setLabels(ax, "kx ($(scaleLabel)⁻¹)", "kx ($(scaleLabel)⁻¹)")

    return fig, ax
end

# Scatterplot showing the values of all G-vectors
function plotkXYVectors(simulation::SimulationDefinition; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    derivedParameters = DerivedParameters(simulation)
    kVectorSet = derivedParameters.kVectorSet

    fig, ax = create2Dfigure(title="K-vectors")

    k₀ = getk₀(kVectorSet)

    nPoints = numHarmonics(kVectorSet)
    kᵢx = Vector{Float64}(undef,nPoints)
    kᵢy = Vector{Float64}(undef,nPoints)
    for iHarmonic in 1:nPoints
        kᵢx[iHarmonic] = kVectorSet.kᵢNorm[iHarmonic][X]*k₀ * scale
        kᵢy[iHarmonic] = kVectorSet.kᵢNorm[iHarmonic][Y]*k₀ * scale
    end

    ax.scatter(kᵢx, kᵢy )
    setLabels(ax, "kx ($(scaleLabel)⁻¹)", "kx ($(scaleLabel)⁻¹)")

    return fig, ax

end
