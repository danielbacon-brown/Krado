
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
    zeroOrderKvector = kVectorSet.kᵢ[zeroOrderIndex]
    
    # Convert to a plottable scale.  Length should be equal to wavelength. 
    λ₀ = getλ₀(kVectorSet.wavenumber)
    k₀ = getk₀(kVectorSet.wavenumber)
    scaledKvector = zeroOrderKvector/norm(zeroOrderKvector)^2 * λ₀ * k₀

    arrow(latticeCenter[X]/scale, latticeCenter[Y]/scale, scaledKvector[X]/scale, scaledKvector[Y]/scale)
    
    return nothing
end


# Adds a quiver for the zero-order k-vectorSet.  Plots at center of lattice.  Length of vector is equal to wavelength.
# TODO: change it to scaled by n?
# function addInjectedKVectorsToPlot(boundaryConditions::InputByOrderBoundaryConditions, layerStack::Vector{<:LayerDefinition}, matCol::MaterialCollection, kVectorSet::KVectorSet, harmonicsSet::HarmonicsSet, lattice::Lattice; scale=μm)
# 
#     latticeCenter = getLatticeCenterCoordinates(lattice)
#     wavenumber = kVectorSet.wavenumber
# 
#     # Convert to a plottable scale.  Length should be equal to wavelength. 
#     λ₀ = getλ₀(wavenumber)
#     k₀ = getk₀(wavenumber)
# 
#     # Get the k-vectors and polarization vectors
#     k3VectorsBottom, k3VectorsTop, P3VectorsBottom, P3VectorsTop = getInjectedKandPvectors( boundaryConditions, harmonicsSet, kVectorSet, layerStack, matCol, wavenumber )
# 
#     for (ϖ,k) in k3VectorsBottom
#         scaledKvector = k/norm(k)^2 * λ₀ * k₀
#         scaledKvectorReal = real(scaledKvector)
#         arrow(latticeCenter[X]/scale, latticeCenter[Y]/scale, scaledKvectorReal[X]/scale, scaledKvectorReal[Y]/scale)
#     end
#     for (ϖ,k) in k3VectorsTop
#         scaledKvector = k/norm(k)^2 * λ₀ * k₀
#         scaledKvectorReal = real(scaledKvector)
#         arrow(latticeCenter[X]/scale, latticeCenter[Y]/scale, scaledKvectorReal[X]/scale, scaledKvectorReal[Y]/scale)
#     end
# 
#     return nothing
# end


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

# plots the 2d components of the injected k-vector
# function plot2DinjectedKVectors(simulationDefinition::SimulationDefinition; scale=μm)
# 
#     derivedParameters = DerivedParameters(simulationDefinition)
# 
# 
#     fig = PyPlot.figure("Plot injected k-vectors", figsize=(5,5))
#     ax = PyPlot.axes()
# 
#     # Plot lattice unit cell
#     addLatticeToPlot(simulationDefinition.lattice; scale=scale)
# 
#     addInjectedKVectorsToPlot(derivedParameters.boundaryConditions, simulationDefinition.layerStack, simulationDefinition.materialCollection, derivedParameters.kVectorSet, derivedParameters.harmonicsSet, simulationDefinition.lattice; scale=μm)
# 
#     setPlotLimitsAroundLattice(simulationDefinition.lattice, ax; scale=scale)
# 
# end



# Adds a quiver for the zero-order k-vectorSet.  Plots at center of lattice.  Length of vector is equal to wavelength.
# function add3DinjectedKVectorsToPlot(simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm)
# 
#     latticeCenter = getLatticeCenterCoordinates(simulationDef.lattice)
#     wavenumber = getWavenumber(derivedParams.boundaryConditions)
#     Zlimits = getLayerStackPlotLimits(simulationDef.layerStack)
# 
#     # Convert to a plottable scale.  Length should be equal to wavelength. 
#     λ₀ = getλ₀(wavenumber)
#     k₀ = getk₀(wavenumber)
# 
#     # Get the k-vectors and polarization vectors
#     k3VectorsBottom, k3VectorsTop, P3VectorsBottom, P3VectorsTop = getInjectedKandPvectors( derivedParams.boundaryConditions, derivedParams.harmonicsSet, derivedParams.kVectorSet, simulationDef.layerStack, simulationDef.materialCollection, wavenumber )
# 
#     for (ϖ,k) in k3VectorsBottom
#         kVectorOrigin = [latticeCenter[X], latticeCenter[Y], Zlimits[1]]
#         addKVectorToPlot(k, kVectorOrigin, wavenumber; scale=μm)
#         # scaledKvector = k/norm(k)^2 * λ₀ * k₀
#         # scaledKvectorReal = real(scaledKvector)
#         # quiver(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, kVectorOrigin[Z]/scale, scaledKvectorReal[X]/scale, scaledKvectorReal[Y]/scale, scaledKvectorReal[Z]/scale)
#     end
#     for (ϖ,k) in k3VectorsTop
#         kVectorOrigin = [latticeCenter[X], latticeCenter[Y], Zlimits[2]]
#         addKVectorToPlot(k, kVectorOrigin, wavenumber; scale=μm)
#         # scaledKvector = k/norm(k)^2 * λ₀ * k₀
#         # scaledKvectorReal = real(scaledKvector)
#         # quiver(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, kVectorOrigin[Z]/scale, scaledKvectorReal[X]/scale, scaledKvectorReal[Y]/scale, scaledKvectorReal[Z]/scale)
#     end
# 
#     return nothing
# end

# function plot3DinjectedKVectors(simulationDef::SimulationDefinition; scale=μm) 
# 
#     derivedParameters = DerivedParameters(simulationDef)
# 
#     fig = figure("3D Injected K-vectors", figsize=(5,5))
#     ax = Axes3D(fig)
# 
#     totalThickness = calcTotalThickness(simulationDef.layerStack)
# 
#     xLimits, yLimits = getLatticePlotLimits(simulationDef.lattice)
#     zLimits = getLayerStackPlotLimits(simulationDef.layerStack)
# 
#     # Plot 3D lattice:
#     plot3Dlattice(simulationDef.lattice, simulationDef.layerStack; scale=μm)
# 
#     add3DinjectedKVectorsToPlot(simulationDef, derivedParameters; scale=μm)
# 
#     setCubicAxes(ax, xLimits, yLimits, zLimits; scale=scale)
# 
# end

function scaleKVectorToWavelength(k::TU3VectorComplex, wavenumber::Wavenumber)
    λ₀ = getλ₀(wavenumber)
    k₀ = getk₀(wavenumber)
    scaledKvector = k/norm(k)^2 * λ₀ * k₀
    return scaledKvector
end

# TODO: is this used?
function addKVectorToPlot(k::TU3VectorComplex, kVectorOrigin::TU3VectorReal, wavenumber::Wavenumber; scale=μm)
    # Scale k-vector so that the length is equal to the wavelength
    scaledKvector = scaleKVectorToWavelength(k, wavenumber)
    scaledKvectorReal = real(scaledKvector)
    scaledKvectorImag = imag(scaledKvector)
    @show k
    @show scaledKvector
    
    # Add to current figure
    quiver(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, kVectorOrigin[Z]/scale, scaledKvectorReal[X]/scale, scaledKvectorReal[Y]/scale, scaledKvectorReal[Z]/scale, color=PVECTORCOLOR, linestyle=PVECTORREALLINESTYLE)
    quiver(kVectorOrigin[X]/scale, kVectorOrigin[Y]/scale, kVectorOrigin[Z]/scale, scaledKvectorImag[X]/scale, scaledKvectorImag[Y]/scale, scaledKvectorImag[Z]/scale, color=PVECTORCOLOR, linestyle=PVECTORIMAGLINESTYLE)
end

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
# function add3DKandPVectorsToPlot( fieldSetXYZ::FieldSetXYZ, injectedOrderIndices::Vector{<:Integer}, simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm, Escale = 1)
    
    latticeCenter = getLatticeCenterCoordinates(simulationDef.lattice)
    Zlimits = getLayerStackPlotLimits(simulationDef.layerStack)
    wavenumber = getWavenumber(simulationDef)
    # nHarmonics = numHarmonics(derivedParams)
    
    if side == BOTTOM # bottom layer       
        n = getn( first(simulationDef.layerStack), simulationDef.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[BOTTOMINDEX]] )
    else # top layer
        n = getn( last(simulationDef.layerStack), simulationDef.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[TOPINDEX]] )
    end
    # if fieldSetXYZ.isForward # bottom layer       
    #     n = getn( first(simulationDef.layerStack), simulationDef.materialCollection, wavenumber)
    #     kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[BOTTOMINDEX]] )
    # else # top layer
    #     n = getn( last(simulationDef.layerStack), simulationDef.materialCollection, wavenumber)
    #     kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[TOPINDEX]] )
    # end
    
    # injectedOrderIndices = getInjectedOrderIndices(fieldSetXYZ)
    
    # add3DKandPVectorsToPlot(fieldSetXYZ, injectedOrderIndices, kVectorOrigin, wavenumber; scale=scale, Escale=Escale)
    
    for ϖindex in injectedOrderIndices
            kXY = getkXY(derivedParams.kVectorSet, ϖindex)
            kXYZ = kXYtokXYZ(kXY, n, wavenumber, fieldSetXYZ.isForward)
            E = fieldSetXYZ.fields[ϖindex,:]
            add3DKandPVectorToPlot(kXYZ, E, kVectorOrigin, wavenumber; scale=scale, Escale=Escale)
    end
    # for ϖindex in 1:nHarmonics
    #     # If the field is non-zero
    #     if sum(abs.(fieldSetXYZ.fields[ϖindex,:])) > 0
    #         kXY = getkXY(derivedParams.kVectorSet, ϖindex)
    #         kXYZ = kXYtokXYZ(kXY, n, wavenumber, fieldSetXYZ.isForward)
    #         E = fieldSetXYZ.fields[ϖindex,:]
    #         add3DKandPVectorsToPlot(kXYZ, E, kVectorOrigin, wavenumber; scale=scale, Escale=Escale)
    #     end
    # end

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
    # nHarmonics = numHarmonics(derivedParams)
    
    if fieldSetXYZ.isForward # bottom layer       
        n = getn( first(simulationDef.layerStack), simulationDef.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[BOTTOMINDEX]] )
    else # top layer
        n = getn( last(simulationDef.layerStack), simulationDef.materialCollection, wavenumber)
        kVectorOrigin = _3VectorFloat([latticeCenter[X], latticeCenter[Y], Zlimits[TOPINDEX]] )
    end
    
    for ϖindex in injectedOrderIndices
        # If the field is non-zero
        # if sum(abs.(fieldSetXYZ.fields[ϖindex,:])) > 0
        kXY = getkXY(derivedParams.kVectorSet, ϖindex)
        kXYZ = kXYtokXYZ(kXY, n, wavenumber, fieldSetXYZ.isForward)
        E = fieldSetXYZ.fields[ϖindex,:]
        add2DKandPVectorToPlot(kXYZ, E, kVectorOrigin, wavenumber; scale=scale, Escale=Escale)
        # end
    end
    # for ϖindex in 1:nHarmonics
    #     # If the field is non-zero
    #     if sum(abs.(fieldSetXYZ.fields[ϖindex,:])) > 0
    #         kXY = getkXY(derivedParams.kVectorSet, ϖindex)
    #         kXYZ = kXYtokXYZ(kXY, n, wavenumber, fieldSetXYZ.isForward)
    #         E = fieldSetXYZ.fields[ϖindex,:]
    #         add2DKandPVectorsToPlot(kXYZ, E, kVectorOrigin, wavenumber; scale=scale, Escale=Escale)
    #     end
    # end

end


# Plots the K and P vectors for injected modes.  Plots at center of lattice.  Length of vector is equal to wavelength.
function add3DinjectedKandPVectorsToPlot(simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm, Escale=1)

    # latticeCenter = getLatticeCenterCoordinates(simulationDef.lattice)
    wavenumber = getWavenumber(derivedParams.boundaryConditions)
    # Zlimits = getLayerStackPlotLimits(simulationDef.layerStack)
    # nHarmonics = numHarmonics(derivedParams)
    
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


# TODO: Implement!.  Needs the output modes as well as input.
# Plots the K and P vectors for harmonic modes in the list.  Plots at center of lattice.  Length of vector is equal to wavelength.
function add3DlistedKandPVectorsToPlot(inputFields::InputFields, outputFields::OutputFields, bottomOrders::Vector{_2VectorInt}, topOrders::Vector{_2VectorInt}, simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm, Escale=1)
    # println("here 1")
    # latticeCenter = getLatticeCenterCoordinates(simulationDef.lattice)
    wavenumber = getWavenumber(derivedParams.boundaryConditions)
    # Zlimits = getLayerStackPlotLimits(simulationDef.layerStack)
    # nHarmonics = numHarmonics(derivedParams)
    
    # Create fields in terms of SP
    # fieldSetSPbottom, fieldSetSPtop = calcSPinputFields( derivedParams.boundaryConditions, derivedParams.harmonicsSet)
    # Convert SP fields to xyz
    # fieldSetXYZbottom, fieldSetXYZtop = convertSPinputFieldsToXYZ( fieldSetSPbottom, fieldSetSPtop, derivedParams.kVectorSet, simulationDef.layerStack, simulationDef.materialCollection, wavenumber)
    
    # Convert inputFields to XYZ field sets
    # inputFieldsBottom = inputFields.bottom
    inputFieldsBottomXYZ = convertFieldSetStackToFieldSetXYZ(inputFields.bottom, derivedParams.kVectorSet, derivedParams.kzBottom)
    inputFieldsTopXYZ = convertFieldSetStackToFieldSetXYZ(inputFields.top, derivedParams.kVectorSet, derivedParams.kzTop)
    outputFieldsBottomXYZ = convertFieldSetStackToFieldSetXYZ(outputFields.bottom, derivedParams.kVectorSet, derivedParams.kzBottom)
    outputFieldsTopXYZ = convertFieldSetStackToFieldSetXYZ(outputFields.top, derivedParams.kVectorSet, derivedParams.kzTop)
    
    
    
    
    # for ϖ in bottomOrders
        # index = getOrderIndex(derivedParams.harmonicsSet, ϖ)
        
    # end
    # bottom input
    bottomIndices = [ getOrderIndex(derivedParams.harmonicsSet, ϖ) for ϖ in bottomOrders]
    add3DKandPVectorsToPlot( inputFieldsBottomXYZ, bottomIndices, BOTTOM, simulationDef, derivedParams; scale=scale, Escale=Escale)
    # top input
    topIndices = [ getOrderIndex(derivedParams.harmonicsSet, ϖ) for ϖ in topOrders]
    add3DKandPVectorsToPlot( inputFieldsTopXYZ, topIndices, TOP, simulationDef, derivedParams; scale=scale, Escale=Escale)

    # bottom output
    bottomIndices = [ getOrderIndex(derivedParams.harmonicsSet, ϖ) for ϖ in bottomOrders]
    add3DKandPVectorsToPlot( outputFieldsBottomXYZ, bottomIndices, BOTTOM, simulationDef, derivedParams; scale=scale, Escale=Escale)
    # @show outputFieldsBottomXYZ.isForward
    # # top input
    topIndices = [ getOrderIndex(derivedParams.harmonicsSet, ϖ) for ϖ in topOrders]
    add3DKandPVectorsToPlot( outputFieldsTopXYZ, topIndices, TOP, simulationDef, derivedParams; scale=scale, Escale=Escale)


    
    # injectedOrderIndicesBottom = getInjectedOrderIndices(fieldSetXYZbottom)
    # add3DKandPVectorsToPlot( fieldSetXYZbottom, injectedOrderIndicesBottom, simulationDef, derivedParams; scale=scale, Escale=Escale)
    # injectedOrderIndicesTop = getInjectedOrderIndices(fieldSetXYZtop)
    # add3DKandPVectorsToPlot( fieldSetXYZtop, injectedOrderIndicesTop, simulationDef, derivedParams; scale=scale, Escale=Escale)
        
    return nothing
end
# add3DlistedKandPVectorsToPlot(inputFields::InputFields, outputFields::OutputFields, bottomOrders, topOrders, simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm, Escale=1) = add3DlistedKandPVectorsToPlot(inputFields, outputFields, convert(Vector{_2VectorInt},bottomOrders), convert(Vector{_2VectorInt},topOrders), simulationDef, derivedParams; scale=scale, Escale=Escale)
# add3DlistedKandPVectorsToPlot(inputFields::InputFields, outputFields::OutputFields, bottomOrders::Vector{T}, topOrders::Vector{W}, simulationDef::SimulationDefinition, derivedParams::DerivedParameters; scale=μm, Escale=1) where {T<:Union{_2VectorInt, Vector{<:Integer}}, W<:Union{_2VectorInt, Vector{<:Integer}}}    = add3DlistedKandPVectorsToPlot(inputFields, outputFields, convert(Vector{_2VectorInt},bottomOrders), convert(Vector{_2VectorInt},topOrders), simulationDef, derivedParams; scale=scale, Escale=Escale)

function plot3DinjectedKandPVectors(simulationDef::SimulationDefinition; scale=μm, Escale=1) 
    
    derivedParameters = DerivedParameters(simulationDef)
    

    fig = figure("3D Injected Modes", figsize=(5,5))
    ax = Axes3D(fig)

    totalThickness = calcTotalThickness(simulationDef.layerStack)
    
    xLimits, yLimits = getLatticePlotLimits(simulationDef.lattice)
    zLimits = getLayerStackPlotLimits(simulationDef.layerStack)
    
    # Plot 3D lattice:
    plot3Dlattice(simulationDef.lattice, simulationDef.layerStack; scale=scale)

    add3DinjectedKandPVectorsToPlot(simulationDef, derivedParameters; scale=scale, Escale=Escale)
        
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

    # latticeCenter = getLatticeCenterCoordinates(simulationDef.lattice)
    wavenumber = getWavenumber(derivedParams.boundaryConditions)
    # Zlimits = getLayerStackPlotLimits(simulationDef.layerStack)
    nHarmonics = numHarmonics(derivedParams)
    
    # Create fields in terms of SP
    fieldSetSPbottom, fieldSetSPtop = calcSPinputFields( derivedParams.boundaryConditions, derivedParams.harmonicsSet)
    # Convert SP fields to xyz
    fieldSetXYZbottom, fieldSetXYZtop = convertSPinputFieldsToXYZ( fieldSetSPbottom, fieldSetSPtop, derivedParams.kVectorSet, simulationDef.layerStack, simulationDef.materialCollection, wavenumber)
    
    injectedOrderIndicesBottom = getInjectedOrderIndices(fieldSetXYZbottom)
    add2DKandPVectorsToPlot( fieldSetXYZbottom, injectedOrderIndicesBottom, simulationDef, derivedParams; scale=scale, Escale=Escale)
    
    injectedOrderIndicesTop = getInjectedOrderIndices(fieldSetXYZtop)
    add2DKandPVectorsToPlot( fieldSetXYZtop, injectedOrderIndicesTop, simulationDef, derivedParams; scale=scale, Escale=Escale)
        
        # arrow(latticeCenter[X]/scale, latticeCenter[Y]/scale, scaledKvectorReal[X]/scale, scaledKvectorReal[Y]/scale)
        
    return nothing
end
add2DinjectedKandPVectorsToPlot(simulationDef::SimulationDefinition; scale=μm, Escale=1) = add3DinjectedKandPVectorsToPlot(simulationDef, DerivedParameters(simulationDef); scale=scale, Escale=Escale)

function plot2DinjectedKandPVectors(simulationDef::SimulationDefinition; scale=μm, Escale=1) 
    
    derivedParameters = DerivedParameters(simulationDef)
    

    fig = figure("2D Injected Modes", figsize=(5,5))
    ax = PyPlot.axes()
    # ax = Axes2D(fig)
    
    addLatticeToPlot(simulationDef.lattice; scale=scale)
        
    # totalThickness = calcTotalThickness(simulationDef.layerStack)
    
    xLimits, yLimits = getLatticePlotLimits(simulationDef.lattice)
    zLimits = getLayerStackPlotLimits(simulationDef.layerStack)
    
    # Plot 3D lattice:
    # plot3Dlattice(simulationDef.lattice, simulationDef.layerStack; scale=scale)

    add2DinjectedKandPVectorsToPlot(simulationDef, derivedParameters; scale=scale, Escale=Escale)
        
    # setCubicAxes(ax, xLimits, yLimits, zLimits; scale=scale)
    setPlotLimitsAroundLattice(simulationDef.lattice, ax; scale=scale)  

end
