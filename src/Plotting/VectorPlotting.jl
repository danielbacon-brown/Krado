
function add3DkVectorToPlot(ax, kVectorOrigin::TU3VectorComplex, scaledKvector::TU3VectorComplex; scale=1)
    scaledKvectorReal = real(scaledKvector)
    scaledKvectorImag = imag(scaledKvector)
    ax.quiver(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, kVectorOrigin[Z]/scale, scaledKvectorReal[X]/scale, scaledKvectorReal[Y]/scale, scaledKvectorReal[Z]/scale, color=KVECTORCOLOR, linestyle=KVECTORREALLINESTYLE)
    ax.quiver(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, kVectorOrigin[Z]/scale, scaledKvectorImag[X]/scale, scaledKvectorImag[Y]/scale, scaledKvectorImag[Z]/scale, color=KVECTORCOLOR, linestyle=KVECTORIMAGLINESTYLE)
end
function add3DpolarizationVectorToPlot(ax, Porigin::TU3VectorComplex, P::TU3VectorComplex; scale=1, Escale=1)

    PReal = real(P)
    PImag = imag(P)

    ax.quiver(Porigin[X]/scale, Porigin[Y]/scale, Porigin[Z]/scale, PReal[X]*Escale, PReal[Y]*Escale, PReal[Z]*Escale, color=PVECTORCOLOR, linestyle=PVECTORREALLINESTYLE)
    ax.quiver(Porigin[X]/scale, Porigin[Y]/scale, Porigin[Z]/scale, PImag[X]*Escale, PImag[Y]*Escale, PImag[Z]*Escale, color=PVECTORCOLOR, linestyle=PVECTORIMAGLINESTYLE)
end


function add2DkVectorToPlot(ax, kVectorOrigin::TU3VectorComplex, k::TU3VectorComplex; scale=1)
    ax.arrow(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, k[X]/scale, k[Y]/scale, color=KVECTORCOLOR, linestyle=KVECTORREALLINESTYLE, width=KVECTORWIDTH/scale, length_includes_head=true)
end

function add2DpolarizationVectorToPlot(ax, Porigin::TU3VectorComplex, P::TU3VectorComplex; scale=1, Escale=1)
    PReal = real(P)
    PImag = imag(P)
    ax.arrow(Porigin[X]/scale, Porigin[Y]/scale, PReal[X]*Escale, PReal[Y]*Escale, color=PVECTORCOLOR, linestyle=KVECTORREALLINESTYLE, width=PVECTORWIDTH/scale, length_includes_head=true)
    ax.arrow(Porigin[X]/scale, Porigin[Y]/scale, PImag[X]*Escale, PImag[Y]*Escale, color=PVECTORCOLOR, linestyle=PVECTORIMAGLINESTYLE, width=PVECTORWIDTH/scale, length_includes_head=true)
end

# Adds a quiver for the zero-order k-vectorSet.  Plots at center of lattice.  Length of vector is equal to wavelength.
# TODO: change it to scaled by n?
# TODO: REPLACE WITH VECTOR PLOTTING BY ORDER?
function addZeroOrderKVectorToPlot(ax, kVectorSet::KVectorSet, harmonicsSet::HarmonicsSet, lattice::Lattice; scale=1)

    latticeCenter = getLatticeCenterCoordinates(lattice)

    zeroOrderIndex = getOrderIndex(harmonicsSet,_2VectorInt(0,0))
    zeroOrderKvector = kVectorSet.kᵢNorm[zeroOrderIndex] * getk₀(kVectorSet.wavenumber)

    # Convert to a plottable scale.  Length should be equal to wavelength.
    λ₀ = getλ₀(kVectorSet.wavenumber)
    k₀ = getk₀(kVectorSet.wavenumber)
    scaledKvector = zeroOrderKvector/norm(zeroOrderKvector)^2 * λ₀ * k₀

    ax.arrow(latticeCenter[X]/scale, latticeCenter[Y]/scale, scaledKvector[X]/scale, scaledKvector[Y]/scale)

    return nothing
end


# plots the 2d components of the 0-order k-vector
# TODO REPLaCE with plot by order
function plot2DZeroOrderKVector(simulationDefinition::SimulationDefinition; scale=1)

    derivedParameters = DerivedParameters(simulationDefinition)

    fig, ax = create2Dfigure(title="Plot 0-order k-vector")

    addLatticeToPlot(ax, simulationDefinition.lattice; scale=scale)

    addZeroOrderKVectorToPlot(ax, derivedParameters.kVectorSet, derivedParameters.harmonicsSet, simulationDefinition.lattice; scale=scale)

    setPlotLimitsAroundLattice(ax, simulationDefinition.lattice; scale=scale)

    return fig, ax
end


function scaleKVectorToWavelength(k::TU3VectorComplex, wavenumber::Wavenumber)
    λ₀ = getλ₀(wavenumber)
    k₀ = getk₀(wavenumber)
    scaledKvector = k/norm(k)^2 * λ₀ * k₀
    return scaledKvector
end

function add3DKandPVectorToPlot(ax, k::TU3VectorComplex, P::TU3VectorComplex, kVectorOrigin::TU3VectorReal, wavenumber::Wavenumber; scale=1, Escale=1)

    # Scale k-vector so that the length is equal to the wavelength
    scaledKvector = scaleKVectorToWavelength(k, wavenumber)
    scaledKvectorReal = real(scaledKvector)
    Porigin = kVectorOrigin + scaledKvectorReal

    # Add to current figure
    add3DkVectorToPlot(ax, kVectorOrigin, scaledKvector; scale=scale)

    add3DpolarizationVectorToPlot(ax, Porigin, P; scale=scale, Escale=Escale)

end


function add3DKandPVectorsToPlot( ax, fieldSetXYZ::FieldSetXYZ, injectedOrderIndices::Vector{<:Integer}, direction::Bool, side::Bool, simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters; scale=1, Escale=1)

    latticeCenter = getLatticeCenterCoordinates(simulationDefinition.lattice)
    Zlimits = getLayerStackPlotLimits(simulationDefinition.layerStack)
    wavenumber = getWavenumber(simulationDefinition)

    if side == BOTTOM # bottom layer
        n = getn( first(simulationDefinition.layerStack), simulationDefinition.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[BOTTOMINDEX]] )
    else # top layer
        n = getn( last(simulationDefinition.layerStack), simulationDefinition.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[TOPINDEX]] )
    end

    for ϖindex in injectedOrderIndices
            kXY = getkXYnorm(derivedParameters.kVectorSet, ϖindex) * getk₀(derivedParameters.kVectorSet)
            kXYZ = kXYtokXYZ(kXY, n, wavenumber, fieldSetXYZ.isForward)
            E = fieldSetXYZ.fields[ϖindex,:]
            add3DKandPVectorToPlot(ax, kXYZ, E, kVectorOrigin, wavenumber; scale=scale, Escale=Escale)
    end

end


function getInjectedOrderIndices(fieldSetXYZ::FieldSetXYZ)
    nHarmonics = numHarmonics(fieldSetXYZ)
    injectedOrderIndices = Int64[]
    for ϖindex in 1:nHarmonics
        if sum(abs.(fieldSetXYZ.fields[ϖindex,:])) > 0 # If the field is non-zero
            push!(injectedOrderIndices, ϖindex)
        end
    end
    return injectedOrderIndices
end


function add2DKandPVectorsToPlot( ax, fieldSetXYZ::FieldSetXYZ, injectedOrderIndices::Vector{<:Integer}, simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters; scale=1, Escale = 1)

    latticeCenter = getLatticeCenterCoordinates(simulationDefinition.lattice)
    Zlimits = getLayerStackPlotLimits(simulationDefinition.layerStack)
    wavenumber = getWavenumber(simulationDefinition)

    if fieldSetXYZ.isForward # bottom layer
        n = getn( first(simulationDefinition.layerStack), simulationDefinition.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[BOTTOMINDEX]] )
    else # top layer
        n = getn( last(simulationDefinition.layerStack), simulationDefinition.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[TOPINDEX]] )
    end

    for ϖindex in injectedOrderIndices
        # If the field is non-zero
        kXY = getkXYnorm(derivedParameters.kVectorSet, ϖindex) * getk₀(derivedParameters.kVectorSet)
        kXYZ = kXYtokXYZ(kXY, n, wavenumber, fieldSetXYZ.isForward)
        E = fieldSetXYZ.fields[ϖindex,:]
        add2DKandPVectorToPlot(ax, kXYZ, E, kVectorOrigin, wavenumber; scale=scale, Escale=Escale)
    end
end

# TODO: make this calculation of SP input fields and inputFieldStack the same thing.
# Plots the K and P vectors for injected modes.  Plots at center of lattice.  Length of vector is equal to wavelength.
function add3DinjectedKandPVectorsToPlot(ax, simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters; scale=1, Escale=1)

    wavenumber = getWavenumber(derivedParameters.boundaryConditions)


    # Create fields in terms of SP
    fieldSetSPbottom, fieldSetSPtop = calcSPinputFields( derivedParameters.boundaryConditions, derivedParameters.harmonicsSet)
    # Convert SP fields to xyz
    fieldSetXYZbottom, fieldSetXYZtop = convertSPinputFieldsToXYZ( fieldSetSPbottom, fieldSetSPtop, derivedParameters.kVectorSet, simulationDefinition.layerStack, simulationDefinition.materialCollection, wavenumber)

    injectedOrderIndicesBottom = getInjectedOrderIndices(fieldSetXYZbottom)
    add3DKandPVectorsToPlot( ax, fieldSetXYZbottom, injectedOrderIndicesBottom, FORWARD, BOTTOM, simulationDefinition, derivedParameters; scale=scale, Escale=Escale)
    injectedOrderIndicesTop = getInjectedOrderIndices(fieldSetXYZtop)
    add3DKandPVectorsToPlot( ax, fieldSetXYZtop, injectedOrderIndicesTop, BACKWARD, TOP, simulationDefinition, derivedParameters; scale=scale, Escale=Escale)

    return nothing
end
add3DinjectedKandPVectorsToPlot(ax, simulationDefinition::SimulationDefinition; scale=1, Escale=1) = add3DinjectedKandPVectorsToPlot(ax, simulationDefinition, DerivedParameters(simulationDefinition); scale=scale, Escale=Escale)

# Plots listed KP vectors. Useful for visualizing diffraction amplitudes and polarization.
function add3DlistedKandPVectorsToPlot(ax, inputFields, outputFields, bottomOrders, topOrders, simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters; scale=1, Escale=1)
    kVectorSet = derivedParameters.kVectorSet
    harmonicsSet = derivedParameters.harmonicsSet
    wavenumber = getWavenumber(derivedParameters.boundaryConditions)

    bottomOrderIndices = [ getOrderIndex(harmonicsSet, ϖ) for ϖ in bottomOrders]
    topOrderIndices = [ getOrderIndex(harmonicsSet, ϖ) for ϖ in topOrders]

    # Create fields in terms of SP
    inputFieldSetXYZbottom = convertFieldSetStackToXYZ(inputFields.bottom, kVectorSet, -derivedParameters.kzNormBottom)
    inputFieldSetXYZtop = convertFieldSetStackToXYZ(inputFields.top, kVectorSet, derivedParameters.kzNormTop)
    outputFieldSetXYZbottom = convertFieldSetStackToXYZ(outputFields.bottom, kVectorSet, derivedParameters.kzNormBottom)
    outputFieldSetXYZtop = convertFieldSetStackToXYZ(outputFields.top, kVectorSet, derivedParameters.kzNormTop)

    # Convert SP fields to xyz
    # INPUT
    add3DKandPVectorsToPlot( ax, inputFieldSetXYZbottom, bottomOrderIndices, FORWARD, BOTTOM, simulationDefinition, derivedParameters; scale=scale, Escale=Escale)
    add3DKandPVectorsToPlot( ax, inputFieldSetXYZtop, topOrderIndices, BACKWARD, TOP, simulationDefinition, derivedParameters; scale=scale, Escale=Escale)
    # OUTPUT
    add3DKandPVectorsToPlot( ax, outputFieldSetXYZbottom, bottomOrderIndices, BACKWARD, BOTTOM, simulationDefinition, derivedParameters; scale=scale, Escale=Escale)
    add3DKandPVectorsToPlot( ax, outputFieldSetXYZtop, topOrderIndices, FORWARD, TOP, simulationDefinition, derivedParameters; scale=scale, Escale=Escale)

end
add3DlistedKandPVectorsToPlot(ax, inputFields, outputFields, bottomOrders, topOrders, simulationDefinition::SimulationDefinition; scale=1, Escale=1) = add3DlistedKandPVectorsToPlot(ax, inputFields, outputFields, bottomOrders, topOrders, simulationDefinition, DerivedParameters(simulationDefinition); scale=1, Escale=1)




function plot3DinjectedKandPVectors(simulationDefinition::SimulationDefinition; scale=1, Escale=1)

    derivedParameters = DerivedParameters(simulationDefinition)

    # Plot 3D lattice:
    fig, ax = plot3Dlattice(simulationDefinition.lattice, simulationDefinition.layerStack; scale=scale, title = "3D Injected Modes")

    add3DinjectedKandPVectorsToPlot(ax, simulationDefinition, derivedParameters; scale=scale, Escale=Escale)

    set3DplotLimits(ax, simulationDefinition.lattice, simulationDefinition.layerStack; scale=scale)

    return fig, ax
end




function add2DKandPVectorToPlot(ax, k::TU3VectorComplex, P::TU3VectorComplex, kVectorOrigin::TU3VectorReal, wavenumber::Wavenumber; scale=1, Escale=1)

    # Scale k-vector so that the length is equal to the wavelength
    scaledKvector = scaleKVectorToWavelength(k, wavenumber)
    scaledKvectorReal = real(scaledKvector)
    scaledKvectorImag = imag(scaledKvector)
    Porigin = kVectorOrigin + scaledKvectorReal

    # Add to current figure
    add2DkVectorToPlot(ax, kVectorOrigin, scaledKvectorReal; scale=scale)
    add2DpolarizationVectorToPlot(ax, Porigin, P; scale=scale, Escale=Escale)
end

# Adds a quiver for the zero-order k-vectorSet.  Plots at center of lattice.  Length of vector is equal to wavelength.
function add2DinjectedKandPVectorsToPlot(ax, simulationDefinition::SimulationDefinition, derivedParameters::DerivedParameters; scale=1, Escale=1)

    wavenumber = getWavenumber(derivedParameters.boundaryConditions)
    nHarmonics = numHarmonics(derivedParameters)

    # Create fields in terms of SP
    fieldSetSPbottom, fieldSetSPtop = calcSPinputFields( derivedParameters.boundaryConditions, derivedParameters.harmonicsSet)
    # Convert SP fields to xyz
    fieldSetXYZbottom, fieldSetXYZtop = convertSPinputFieldsToXYZ( fieldSetSPbottom, fieldSetSPtop, derivedParameters.kVectorSet, simulationDefinition.layerStack, simulationDefinition.materialCollection, wavenumber)

    injectedOrderIndicesBottom = getInjectedOrderIndices(fieldSetXYZbottom)
    add2DKandPVectorsToPlot( ax, fieldSetXYZbottom, injectedOrderIndicesBottom, simulationDefinition, derivedParameters; scale=scale, Escale=Escale)

    injectedOrderIndicesTop = getInjectedOrderIndices(fieldSetXYZtop)
    add2DKandPVectorsToPlot( ax, fieldSetXYZtop, injectedOrderIndicesTop, simulationDefinition, derivedParameters; scale=scale, Escale=Escale)

end
add2DinjectedKandPVectorsToPlot(ax, simulationDefinition::SimulationDefinition; scale=1, Escale=1) = add3DinjectedKandPVectorsToPlot(ax, simulationDefinition, DerivedParameters(simulationDefinition); scale=scale, Escale=Escale)

# CURRENTLY NOT USED
function plot2DinjectedKandPVectors(simulationDefinition::SimulationDefinition; scale=1, Escale=1)

    derivedParameters = DerivedParameters(simulationDefinition)

    fig, ax = create2Dfigure(title="2D Injected Modes")

    addLatticeToPlot(ax, simulationDefinition.lattice; scale=scale)

    add2DinjectedKandPVectorsToPlot(ax, simulationDefinition, derivedParameters; scale=scale, Escale=Escale)

    setPlotLimitsAroundLattice(ax, simulationDefinition.lattice; scale=scale)

    return fig, ax
end
