
function add3DkVectorToPlot(kVectorOrigin::TU3VectorComplex, scaledKvector::TU3VectorComplex; scale=μm)
    scaledKvectorReal = real(scaledKvector)
    scaledKvectorImag = imag(scaledKvector)
    quiver(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, kVectorOrigin[Z]/scale, scaledKvectorReal[X]/scale, scaledKvectorReal[Y]/scale, scaledKvectorReal[Z]/scale, color=KVECTORCOLOR, linestyle=KVECTORREALLINESTYLE)
    quiver(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, kVectorOrigin[Z]/scale, scaledKvectorImag[X]/scale, scaledKvectorImag[Y]/scale, scaledKvectorImag[Z]/scale, color=KVECTORCOLOR, linestyle=KVECTORIMAGLINESTYLE)
end
function add3DpolarizationVectorToPlot(Porigin::TU3VectorComplex, P::TU3VectorComplex; scale=μm, Escale=1)

    PReal = real(P)
    PImag = imag(P)

    quiver(Porigin[X]/scale, Porigin[Y]/scale, Porigin[Z]/scale, PReal[X]*Escale, PReal[Y]*Escale, PReal[Z]*Escale, color=PVECTORCOLOR, linestyle=PVECTORREALLINESTYLE)
    quiver(Porigin[X]/scale, Porigin[Y]/scale, Porigin[Z]/scale, PImag[X]*Escale, PImag[Y]*Escale, PImag[Z]*Escale, color=PVECTORCOLOR, linestyle=PVECTORIMAGLINESTYLE)
end


function add2DkVectorToPlot(kVectorOrigin::TU3VectorComplex, k::TU3VectorComplex; scale=μm)
    arrow(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, k[X]/scale, k[Y]/scale, color=KVECTORCOLOR, linestyle=KVECTORREALLINESTYLE, width=KVECTORWIDTH/scale, length_includes_head=true)
end

function add2DpolarizationVectorToPlot(Porigin::TU3VectorComplex, P::TU3VectorComplex; scale=μm, Escale=1)
    PReal = real(P)
    PImag = imag(P)
    arrow(Porigin[X]/scale, Porigin[Y]/scale, PReal[X]*Escale, PReal[Y]*Escale, color=PVECTORCOLOR, linestyle=KVECTORREALLINESTYLE, width=PVECTORWIDTH/scale, length_includes_head=true)
    arrow(Porigin[X]/scale, Porigin[Y]/scale, PImag[X]*Escale, PImag[Y]*Escale, color=PVECTORCOLOR, linestyle=PVECTORIMAGLINESTYLE, width=PVECTORWIDTH/scale, length_includes_head=true)
end

# Adds a quiver for the zero-order k-vectorSet.  Plots at center of lattice.  Length of vector is equal to wavelength.
# TODO: change it to scaled by n?
# TODO: REPLACE WITH VECTOR PLOTTING BY ORDER?
function addZeroOrderKVectorToPlot(kVectorSet::KVectorSet, harmonicsSet::HarmonicsSet, lattice::Lattice; scale=μm)

    latticeCenter = getLatticeCenterCoordinates(lattice)

    zeroOrderIndex = getOrderIndex(harmonicsSet,_2VectorInt(0,0))
    zeroOrderKvector = kVectorSet.kᵢNorm[zeroOrderIndex] * getk₀(kVectorSet.wavenumber)

    # Convert to a plottable scale.  Length should be equal to wavelength.
    λ₀ = getλ₀(kVectorSet.wavenumber)
    k₀ = getk₀(kVectorSet.wavenumber)
    scaledKvector = zeroOrderKvector/norm(zeroOrderKvector)^2 * λ₀ * k₀

    arrow(latticeCenter[X]/scale, latticeCenter[Y]/scale, scaledKvector[X]/scale, scaledKvector[Y]/scale)

    return nothing
end


# plots the 2d components of the 0-order k-vector
# TODO REPLaCE with plot by order
function plot2DZeroOrderKVector(simulationDefinition::SimulationDefinition; scale=μm)

    derivedParameters = DerivedParameters(simulationDefinition)


    fig = PyPlot.figure("Plot 0-order k-vector", figsize=(5,5))
    ax = PyPlot.axes()

    # Plot lattice unit cell
    addLatticeToPlot(simulationDefinition.lattice; scale=scale)

    addZeroOrderKVectorToPlot(derivedParameters.kVectorSet, derivedParameters.harmonicsSet, simulationDefinition.lattice; scale=scale)

    setPlotLimitsAroundLattice(simulationDefinition.lattice, ax; scale=scale)

end


function scaleKVectorToWavelength(k::TU3VectorComplex, wavenumber::Wavenumber)
    λ₀ = getλ₀(wavenumber)
    k₀ = getk₀(wavenumber)
    scaledKvector = k/norm(k)^2 * λ₀ * k₀
    return scaledKvector
end

# TODO: is this used?
# function addKVectorToPlot(k::TU3VectorComplex, kVectorOrigin::TU3VectorReal, wavenumber::Wavenumber; scale=μm)
#     # Scale k-vector so that the length is equal to the wavelength
#     scaledKvector = scaleKVectorToWavelength(k, wavenumber)
#     scaledKvectorReal = real(scaledKvector)
#     scaledKvectorImag = imag(scaledKvector)
#     # @show k
#     # @show scaledKvector
#
#     # Add to current figure
#     quiver(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, kVectorOrigin[Z]/scale, scaledKvectorReal[X]/scale, scaledKvectorReal[Y]/scale, scaledKvectorReal[Z]/scale, color=PVECTORCOLOR, linestyle=PVECTORREALLINESTYLE)
#     quiver(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, kVectorOrigin[Z]/scale, scaledKvectorImag[X]/scale, scaledKvectorImag[Y]/scale, scaledKvectorImag[Z]/scale, color=PVECTORCOLOR, linestyle=PVECTORIMAGLINESTYLE)
# end

function add3DKandPVectorToPlot(k::TU3VectorComplex, P::TU3VectorComplex, kVectorOrigin::TU3VectorReal, wavenumber::Wavenumber; scale=μm, Escale=1)

    # Scale k-vector so that the length is equal to the wavelength
    scaledKvector = scaleKVectorToWavelength(k, wavenumber)
    scaledKvectorReal = real(scaledKvector)
    Porigin = kVectorOrigin + scaledKvectorReal

    # Add to current figure
    add3DkVectorToPlot(kVectorOrigin, scaledKvector; scale=scale)

    add3DpolarizationVectorToPlot(Porigin, P; scale=scale, Escale=Escale)

end


function add3DKandPVectorsToPlot( fieldSetXYZ::FieldSetXYZ, injectedOrderIndices::Vector{<:Integer}, side::Bool, simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm, Escale=1)

    latticeCenter = getLatticeCenterCoordinates(simulationDef.lattice)
    Zlimits = getLayerStackPlotLimits(simulationDef.layerStack)
    wavenumber = getWavenumber(simulationDef)

    if side == BOTTOM # bottom layer
        n = getn( first(simulationDef.layerStack), simulationDef.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[BOTTOMINDEX]] )
    else # top layer
        n = getn( last(simulationDef.layerStack), simulationDef.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[TOPINDEX]] )
    end

    for ϖindex in injectedOrderIndices
            # kXY = getkXY(derivedParams.kVectorSet, ϖindex)
            kXY = getkXYnorm(derivedParams.kVectorSet, ϖindex) * getk₀(derivedParams.kVectorSet)
            kXYZ = kXYtokXYZ(kXY, n, wavenumber, fieldSetXYZ.isForward)
            E = fieldSetXYZ.fields[ϖindex,:]
            add3DKandPVectorToPlot(kXYZ, E, kVectorOrigin, wavenumber; scale=scale, Escale=Escale)
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


function add2DKandPVectorsToPlot( fieldSetXYZ::FieldSetXYZ, injectedOrderIndices::Vector{<:Integer}, simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm, Escale = 1)

    latticeCenter = getLatticeCenterCoordinates(simulationDef.lattice)
    Zlimits = getLayerStackPlotLimits(simulationDef.layerStack)
    wavenumber = getWavenumber(simulationDef)

    if fieldSetXYZ.isForward # bottom layer
        n = getn( first(simulationDef.layerStack), simulationDef.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[BOTTOMINDEX]] )
    else # top layer
        n = getn( last(simulationDef.layerStack), simulationDef.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[TOPINDEX]] )
    end

    for ϖindex in injectedOrderIndices
        # If the field is non-zero
        kXY = getkXYnorm(derivedParams.kVectorSet, ϖindex) * getk₀(derivedParams.kVectorSet)
        kXYZ = kXYtokXYZ(kXY, n, wavenumber, fieldSetXYZ.isForward)
        E = fieldSetXYZ.fields[ϖindex,:]
        add2DKandPVectorToPlot(kXYZ, E, kVectorOrigin, wavenumber; scale=scale, Escale=Escale)
    end
end

# TODO: make this calculation of SP input fields and inputFieldStack the same thing.
# Plots the K and P vectors for injected modes.  Plots at center of lattice.  Length of vector is equal to wavelength.
function add3DinjectedKandPVectorsToPlot(simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm, Escale=1)

    wavenumber = getWavenumber(derivedParams.boundaryConditions)


    # Create fields in terms of SP
    fieldSetSPbottom, fieldSetSPtop = calcSPinputFields( derivedParams.boundaryConditions, derivedParams.harmonicsSet)
    # Convert SP fields to xyz
    fieldSetXYZbottom, fieldSetXYZtop = convertSPinputFieldsToXYZ( fieldSetSPbottom, fieldSetSPtop, derivedParams.kVectorSet, simulationDef.layerStack, simulationDef.materialCollection, wavenumber)

    injectedOrderIndicesBottom = getInjectedOrderIndices(fieldSetXYZbottom)
    add3DKandPVectorsToPlot( fieldSetXYZbottom, injectedOrderIndicesBottom, BOTTOM, simulationDef, derivedParams; scale=scale, Escale=Escale)
    injectedOrderIndicesTop = getInjectedOrderIndices(fieldSetXYZtop)
    add3DKandPVectorsToPlot( fieldSetXYZtop, injectedOrderIndicesTop, TOP, simulationDef, derivedParams; scale=scale, Escale=Escale)

    return nothing
end
add3DinjectedKandPVectorsToPlot(simulationDef::SimulationDefinition; scale=μm, Escale=1) = add3DinjectedKandPVectorsToPlot(simulationDef, DerivedParameters(simulationDef); scale=scale, Escale=Escale)

# Plots listed KP vectors. Useful for visualizing diffraction amplitudes and polarization.
function add3DlistedKandPVectorsToPlot(inputFields, outputFields, bottomOrders, topOrders, simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm, Escale=1)
# allModeData.inputFields, allModeData.outputFields, bottomOrders, topOrders, simulationDefinition, derivedParameters; scale=μm, Escale = 0.3
    kVectorSet = derivedParams.kVectorSet
    harmonicsSet = derivedParams.harmonicsSet
    wavenumber = getWavenumber(derivedParams.boundaryConditions)

    bottomOrderIndices = [ getOrderIndex(harmonicsSet, ϖ) for ϖ in bottomOrders]
    topOrderIndices = [ getOrderIndex(harmonicsSet, ϖ) for ϖ in topOrders]

    # Create fields in terms of SP
    # fieldSetSPbottom, fieldSetSPtop = calcSPinputFields( derivedParams.boundaryConditions, derivedParams.harmonicsSet)


    inputFieldSetXYZbottom = convertFieldSetStackToXYZ(inputFields.bottom, kVectorSet, derivedParams.kzNormBottom)
    inputFieldSetXYZtop = convertFieldSetStackToXYZ(inputFields.top, kVectorSet, derivedParams.kzNormTop)
    outputFieldSetXYZbottom = convertFieldSetStackToXYZ(outputFields.bottom, kVectorSet, derivedParams.kzNormBottom)
    outputFieldSetXYZtop = convertFieldSetStackToXYZ(outputFields.top, kVectorSet, derivedParams.kzNormTop)

    # Convert SP fields to xyz
    # fieldSetXYZbottom, fieldSetXYZtop = convertSPinputFieldsToXYZ( fieldSetSPbottom, fieldSetSPtop, derivedParams.kVectorSet, simulationDef.layerStack, simulationDef.materialCollection, wavenumber)

    # injectedOrderIndicesBottom = getInjectedOrderIndices(fieldSetXYZbottom)
    add3DKandPVectorsToPlot( inputFieldSetXYZbottom, bottomOrderIndices, BOTTOM, simulationDef, derivedParams; scale=scale, Escale=Escale)
    # add3DKandPVectorsToPlot( inputFieldSetXYZtop, topOrderIndices, TOP, simulationDef, derivedParams; scale=scale, Escale=Escale)
    # injectedOrderIndicesTop = getInjectedOrderIndices(fieldSetXYZtop)
    # add3DKandPVectorsToPlot( fieldSetXYZtop, injectedOrderIndicesTop, TOP, simulationDef, derivedParams; scale=scale, Escale=Escale)

    return nothing
end




function plot3DinjectedKandPVectors(simulationDef::SimulationDefinition; scale=μm, Escale=1)

    derivedParameters = DerivedParameters(simulationDef)


    # fig = figure("3D Injected Modes", figsize=(5,5))
    # ax = Axes3D(fig)
    #
    # totalThickness = calcTotalThickness(simulationDef.layerStack)
    #
    # xLimits, yLimits = getLatticePlotLimits(simulationDef.lattice)
    # zLimits = getLayerStackPlotLimits(simulationDef.layerStack)

    # Plot 3D lattice:
    fig, ax = plot3Dlattice(simulationDef.lattice, simulationDef.layerStack; scale=scale, title = "3D Injected Modes")

    add3DinjectedKandPVectorsToPlot(simulationDef, derivedParameters; scale=scale, Escale=Escale)

    xLimits, yLimits = getLatticePlotLimits(simulationDef.lattice)
    zLimits = getLayerStackPlotLimits(simulationDef.layerStack)
    setCubicAxes(ax, xLimits, yLimits, zLimits; scale=scale)

end




function add2DKandPVectorToPlot(k::TU3VectorComplex, P::TU3VectorComplex, kVectorOrigin::TU3VectorReal, wavenumber::Wavenumber; scale=μm, Escale=1)

    # Scale k-vector so that the length is equal to the wavelength
    scaledKvector = scaleKVectorToWavelength(k, wavenumber)
    scaledKvectorReal = real(scaledKvector)
    scaledKvectorImag = imag(scaledKvector)
    Porigin = kVectorOrigin + scaledKvectorReal

    # Add to current figure
    add2DkVectorToPlot(kVectorOrigin, scaledKvectorReal; scale=scale)
    add2DpolarizationVectorToPlot(Porigin, P; scale=scale, Escale=Escale)
end

# Adds a quiver for the zero-order k-vectorSet.  Plots at center of lattice.  Length of vector is equal to wavelength.
function add2DinjectedKandPVectorsToPlot(simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm, Escale=1)

    wavenumber = getWavenumber(derivedParams.boundaryConditions)
    nHarmonics = numHarmonics(derivedParams)

    # Create fields in terms of SP
    fieldSetSPbottom, fieldSetSPtop = calcSPinputFields( derivedParams.boundaryConditions, derivedParams.harmonicsSet)
    # Convert SP fields to xyz
    fieldSetXYZbottom, fieldSetXYZtop = convertSPinputFieldsToXYZ( fieldSetSPbottom, fieldSetSPtop, derivedParams.kVectorSet, simulationDef.layerStack, simulationDef.materialCollection, wavenumber)

    injectedOrderIndicesBottom = getInjectedOrderIndices(fieldSetXYZbottom)
    add2DKandPVectorsToPlot( fieldSetXYZbottom, injectedOrderIndicesBottom, simulationDef, derivedParams; scale=scale, Escale=Escale)

    injectedOrderIndicesTop = getInjectedOrderIndices(fieldSetXYZtop)
    add2DKandPVectorsToPlot( fieldSetXYZtop, injectedOrderIndicesTop, simulationDef, derivedParams; scale=scale, Escale=Escale)

    return nothing
end
add2DinjectedKandPVectorsToPlot(simulationDef::SimulationDefinition; scale=μm, Escale=1) = add3DinjectedKandPVectorsToPlot(simulationDef, DerivedParameters(simulationDef); scale=scale, Escale=Escale)

function plot2DinjectedKandPVectors(simulationDef::SimulationDefinition; scale=μm, Escale=1)

    derivedParameters = DerivedParameters(simulationDef)


    fig = figure("2D Injected Modes", figsize=(5,5))
    ax = PyPlot.axes()

    addLatticeToPlot(simulationDef.lattice; scale=scale)

    xLimits, yLimits = getLatticePlotLimits(simulationDef.lattice)
    zLimits = getLayerStackPlotLimits(simulationDef.layerStack)

    add2DinjectedKandPVectorsToPlot(simulationDef, derivedParameters; scale=scale, Escale=Escale)

    setPlotLimitsAroundLattice(simulationDef.lattice, ax; scale=scale)

end
