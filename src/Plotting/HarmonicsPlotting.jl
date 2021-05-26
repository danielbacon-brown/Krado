# Harmonics related plotting functions



# Scatterplot showing the (m,n) harmonic orders
function plotHarmonicsSet(simulation::SimulationDefinition)

    derivedParameters = DerivedParameters(simulation)
    harmonicsSet = derivedParameters.harmonicsSet

    fig = PyPlot.figure("Harmonics (m,n)", figsize=(5,5))
    ax = PyPlot.axes()

    nPoints = numHarmonics(harmonicsSet)
    mᵢ = Vector{Int64}(undef,nPoints)
    nᵢ = Vector{Int64}(undef,nPoints)
    for iHarmonic in 1:nPoints
        mᵢ[iHarmonic] = harmonicsSet.mnᵢ[iHarmonic][U]
        nᵢ[iHarmonic] = harmonicsSet.mnᵢ[iHarmonic][V]
    end

    PyPlot.scatter(mᵢ, nᵢ )
    PyPlot.xlabel("m")
    PyPlot.ylabel("n")

    return fig, ax
end


# Scatterplot showing the values of all G-vectors
function plotGvectorSet(simulation::SimulationDefinition; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    derivedParameters = DerivedParameters(simulation)
    gVectorSet = derivedParameters.gVectorSet

    fig = PyPlot.figure("G-vectors", figsize=(5,5))
    ax = PyPlot.axes()

    nPoints = numGvectors(gVectorSet)
    Gᵢx = Vector{Float64}(undef,nPoints)
    Gᵢy = Vector{Float64}(undef,nPoints)
    for iHarmonic in 1:nPoints
        Gᵢx[iHarmonic] = gVectorSet.Gᵢ[iHarmonic][X] * scale
        Gᵢy[iHarmonic] = gVectorSet.Gᵢ[iHarmonic][Y] * scale
    end

    PyPlot.scatter(Gᵢx, Gᵢy )
    PyPlot.xlabel("kx ($(scaleLabel)⁻¹)")
    PyPlot.ylabel("ky ($(scaleLabel)⁻¹)")

    return fig, ax
end

# Scatterplot showing the values of all G-vectors
function plotkXYVectors(simulation::SimulationDefinition; scale=1)

    scaleLabel = LENGTHLABEL[scale]

    derivedParameters = DerivedParameters(simulation)
    kVectorSet = derivedParameters.kVectorSet

    fig = PyPlot.figure("kXY-vectors", figsize=(5,5))
    ax = PyPlot.axes()

    k₀ = getk₀(kVectorSet)

    nPoints = numHarmonics(kVectorSet)
    kᵢx = Vector{Float64}(undef,nPoints)
    kᵢy = Vector{Float64}(undef,nPoints)
    for iHarmonic in 1:nPoints
        kᵢx[iHarmonic] = kVectorSet.kᵢNorm[iHarmonic][X]*k₀ * scale
        kᵢy[iHarmonic] = kVectorSet.kᵢNorm[iHarmonic][Y]*k₀ * scale
    end

    PyPlot.scatter(kᵢx, kᵢy )
    PyPlot.xlabel("kx ($(scaleLabel)⁻¹)")
    PyPlot.ylabel("ky ($(scaleLabel)⁻¹)")

    return fig, ax

end
