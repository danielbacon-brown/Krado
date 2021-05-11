# A layer definition is a structure that describes how the convolution matrix should be calculated
# All have a thickness and a background material

abstract type LayerDefinition end

# Describes a layer that is uniform in the x,y directions.  Contains no patterns/shapes
mutable struct UniformLayerDefinition <: LayerDefinition
    thickness::Float64
    backgroundMaterialName::String

    function UniformLayerDefinition(thickness::Float64, backgroundMaterialName::String)
        return new(thickness, backgroundMaterialName)
    end
end

mutable struct SemiInfiniteLayerDefinition <: LayerDefinition
    backgroundMaterialName::String

    function SemiInfiniteLayerDefinition( backgroundMaterialName::String)
        return new(backgroundMaterialName)
    end
end



# For an FFT calculation of convolution matrix
# Points are positioned in a grid following the lattice periodicity
mutable struct PatternedLayerDefinition <: LayerDefinition
    # Number of points to use in each of the real-space period vectors [Nx, Ny]
    numDivisions::_2VectorInt  # length == 2
    thickness::Float64
    layerPattern::AbstractLayerPattern

    function PatternedLayerDefinition(numDivisions::_2VectorInt, thickness::Real, layerPattern::AbstractLayerPattern)
        return new(convert(_2VectorInt,numDivisions), convert(Float64,thickness), layerPattern)
    end
end
PatternedLayerDefinition(numDivisions, thickness::Real, layerPattern::AbstractLayerPattern) = PatternedLayerDefinition(_2VectorInt(numDivisions), thickness, layerPattern)


# If numDivisions is a scalar, use that as 1D
function PatternedLayerDefinition(numDivisions::Number, thickness::Real, layerPattern::AbstractLayerPattern)
    return PatternedLayerDefinition( _2VectorInt(numDivisions,1), thickness, layerPattern)
end



# TODO: Make more efficient
# Calculate the transform ∫∫x⋅exp(G̅⋅𝐫) for a single G-vector.  Normalizing to the number of coordinates.
function calc∫xexp𝐆𝐫(valueGrid, positionGrid::PositionGridXY, G::_2VectorFloat)

    @assert size(valueGrid) == size(positionGrid.positions)

    gridSize = size(valueGrid)
    gridIndices = getGridIndices(gridSize)

    total = ComplexF64(0)
    for i_u = 1:gridSize[1]
        for i_v = 1:gridSize[2]
            total += valueGrid[i_u, i_v] * exp( -1im * G ⋅ positionGrid.positions[i_u, i_v] )
        end
    end

    return total / (gridSize[1]*gridSize[2])
end

# Calculate the transform ∫∫x⋅exp(G̅⋅𝐫) for all G-vectorDifferences, putting results in a dict with the harmonic mn as the key.
# This corresponds to the 'a' matrix
# valueGrid is a 2D grid of the values at r real-space coordinates.  Does transform using a single G-vector: G.
function calc∫xexpΔ𝐆𝐫Dict(valuesGrid, positionGrid::PositionGridXY, Gvectors::GvectorSet)

    ∫xexpΔ𝐆𝐫 = Dict{_2VectorInt, ComplexF64}()
    for ΔharmonicIndex in 1:numΔGvectors(Gvectors)
        harmonic = Gvectors.harmonicsSet.Δmnᵢⱼ[ΔharmonicIndex]
        ΔG = Gvectors.ΔGᵢⱼ[harmonic]
        push!( ∫xexpΔ𝐆𝐫, harmonic => calc∫xexp𝐆𝐫(valuesGrid, positionGrid, ΔG) )
    end
    return ∫xexpΔ𝐆𝐫
end

# Input: dict describing ∫ϵ⋅exp(Δ𝐆⋅𝐫)
#    and the set describing the harmonics to use
# Output: 2D grid describing the convolution matrix (e.g. ⟦ϵ⟧  or ⟦μ⟧ )
# Iterate over every combination of G-vectors and grab appropriate value from dict
function assembleConvolutionMatrix( valuesByΔϖ::Dict{_2VectorInt,ComplexF64}, harmonicsSet::HarmonicsSet  )

    # TODO: Replace w nHarmonics

    # Refers to either Cϵ, ⟦ϵ⟧  or Cμ, ⟦μ⟧
    Cϵᵢⱼ = Array{ComplexF64,2}(undef, (numHarmonics(harmonicsSet), numHarmonics(harmonicsSet)) )

    for i_G = 1:numHarmonics(harmonicsSet)
        for j_G = 1:numHarmonics(harmonicsSet)
            Cϵᵢⱼ[i_G,j_G] = valuesByΔϖ[ harmonicsSet.mnᵢ[i_G] - harmonicsSet.mnᵢ[j_G]]
        end
    end
    return Cϵᵢⱼ
end

# Calculates the convolution matrix for x (either ϵ or μ) by creating a dictionary of results for all possible pairs of G-vectors, then quickly generating the matrix using the dict
function calcConvolutionMatrix( xGrid, positionGrid::PositionGridXY, Gvectors::GvectorSet)
    ∫xexpΔ𝐆𝐫Dict = calc∫xexpΔ𝐆𝐫Dict(xGrid, positionGrid, Gvectors)
    Cxᵢⱼ = assembleConvolutionMatrix( ∫xexpΔ𝐆𝐫Dict, Gvectors.harmonicsSet )
    return Cxᵢⱼ
end


function getPositionϵμGrids( layerDef::PatternedLayerDefinition, lattice::Lattice, matCol::MaterialCollection, wavenumber::Wavenumber)
    # posGrid = calcUniformGridPositions(lattice, layerDef.numDivisions)
    positionGrid = PositionGridXYleftAligned(lattice, layerDef.numDivisions)
    ϵGrid = getϵAtPosition( layerDef.layerPattern, positionGrid, matCol, wavenumber)
    μGrid = getμAtPosition( layerDef.layerPattern, positionGrid, matCol, wavenumber)
    return positionGrid, ϵGrid, μGrid
end

# Calculates the convolution matrices ⟦ϵ⟧  ⟦μ⟧ a the given layer and lattice with the corresponding harmonics and materials
function calcConvolutionMatrices( layerDef::PatternedLayerDefinition, lattice::Lattice, Gvectors::GvectorSet, matCol::MaterialCollection, wavenumber::Wavenumber )

    positionGrid, ϵGrid, μGrid = getPositionϵμGrids(layerDef, lattice, matCol, wavenumber)

    Cϵᵢⱼ = calcConvolutionMatrix( ϵGrid, positionGrid, Gvectors)
    Cμᵢⱼ = calcConvolutionMatrix( μGrid, positionGrid, Gvectors)

    return Cϵᵢⱼ, Cμᵢⱼ
end

# IS THIS USED? - only in IntegrationTest3
# Calculates the convolution matrices ⟦ϵ⟧  ⟦μ⟧ a the given layer and lattice with the corresponding harmonics and materials.  For an unpatterned layer, the convolution matrix is simply the a diagonal matrix with all components being equal to permittivity
function calcConvolutionMatrices( layerDef::UniformLayerDefinition, lattice::Lattice, Gvectors::GvectorSet, matCol::MaterialCollection, wavenumber::Wavenumber )

    ϵ, μ = calc_ϵμ( getMaterial(matCol,layerDef.backgroundMaterialName), wavenumber)

    # @show ϵ
    # ϵ = conj.(ϵ)  # DOESN'T MATTER FOR INTEGRATION TEST 3
    # @show ϵ

    numHarmonics = numGvectors(Gvectors)
    Cϵᵢⱼ = Matrix{ComplexF64}(ϵ*I, numHarmonics, numHarmonics)
    Cμᵢⱼ = Matrix{ComplexF64}(μ*I, numHarmonics, numHarmonics)

    return Cϵᵢⱼ, Cμᵢⱼ
end




# Returns the position grid Array{_2VectorFloat(X,Y),2}
function calcUniformGridPositions(lattice::Lattice, layerDef::PatternedLayerDefinition)
    return calcUniformGridPositions(lattice, layerDef.numDivisions)
end


#shorthand for gettign the refractive index of a uniform layer
function getn(layer::layerT, matCol::MaterialCollection, wavenumber::Wavenumber) where layerT<:Union{SemiInfiniteLayerDefinition, UniformLayerDefinition}
    return convert_ϵ2n(calc_ϵ( getMaterial(matCol,layer.backgroundMaterialName), wavenumber))
end

function getϵμ(layer::layerT, matCol::MaterialCollection, wavenumber::Wavenumber) where layerT<:Union{SemiInfiniteLayerDefinition, UniformLayerDefinition}
    return calc_ϵμ( getMaterial(matCol,layer.backgroundMaterialName), wavenumber)
end


#Calculates the z-components k-vectors of the modes in the top or bottom layers.  Assumed real components are positive.
# TODO: Figure out why n^2, and ϵ*μ works but not conj(ϵ)*conj(μ)
function calckz(kVectorSet::KVectorSet, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection, wavenumber::Wavenumber)
    n = getn(layer, matCol, wavenumber)
    ϵ,μ = getϵμ( layer, matCol, wavenumber)
    # @assert isapprox(n^2, ϵ*μ, rtol=1e-3)  # should be moved to a test.

    # Lecture 7B appears to take the conjugate of Kz
    # Kz = conj( sqrt(I-Kx^2-Ky^2))
    # This causes some tests to fail
    # Lecture 7B:
    # return ComplexF64[ conj(sqrt(getk₀(kVectorSet)^2*conj(ϵ)*conj(μ) - kᵢ[X]^2 - kᵢ[Y]^2))  for kᵢ in kVectorSet.kᵢ]

    # from edmundsj
    #-conj(sqrt(conj(layer.er*layer.ur)*complexIdentity(kx.shape[0]) - kx @ kx - ky @ ky))
    # based on edmundsj:
    # return ComplexF64[ conj(sqrt(getk₀(kVectorSet)^2*conj(ϵ*μ) - kᵢ[X]^2 - kᵢ[Y]^2))  for kᵢ in kVectorSet.kᵢ]
    # # old:
    return ComplexF64[ sqrt(getk₀(kVectorSet)^2*conj(ϵ)*conj(μ) - kᵢ[X]^2 - kᵢ[Y]^2)  for kᵢ in kVectorSet.kᵢ]

end

# K-NORM
#Calculates the z-components k-vectors of the modes in the top or bottom layers.  Assumed real components are positive.
# TODO: Figure out why n^2, and ϵ*μ works but not conj(ϵ)*conj(μ)
# function calckzGap(kVectorSet::KVectorSet, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection, wavenumber::Wavenumber)
#     ϵ,μ  = 1.0,1.0
#     # return ComplexF64[ conj(sqrt(conj(ϵ*μ)*I - kᵢ[X].^2 - kᵢ[Y].^2))
#     return ComplexF64[ conj(sqrt(conj(ϵ*μ) - kᵢ[X]^2 - kᵢ[Y]^2))  for kᵢ in kVectorSet.kᵢ]
# end

# only difference is that the bottom one has a negative value.
# function calckzBottom(kVectorSet::KVectorSet, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection, wavenumber::Wavenumber)
#     n = getn(layer, matCol, wavenumber)
#     ϵ,μ = getϵμ( layer, matCol, wavenumber)
#     @assert isapprox(n^2, ϵ*μ, rtol=1e-3)
#
#     # Lecture 7B appears to take the conjugate of Kz
#     # Kz = conj( sqrt(I-Kx^2-Ky^2))
#     # This causes some tests to fail.  But seems to work when you also take the conjugate in the lambda calculation and
#     #TODO
#     return -conj( ComplexF64[ sqrt(getk₀(kVectorSet)^2*conj(ϵ)*conj(μ) - kᵢ[X]^2 - kᵢ[Y]^2)  for kᵢ in kVectorSet.kᵢ] )  #TODO
#
#     # old:
#     # return -ComplexF64[ sqrt(getk₀(kVectorSet)^2*conj(ϵ)*conj(μ) - kᵢ[X]^2 - kᵢ[Y]^2)  for kᵢ in kVectorSet.kᵢ]
#
# end
# function calckzTop(kVectorSet::KVectorSet, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection, wavenumber::Wavenumber)
#     n = getn(layer, matCol, wavenumber)
#     ϵ,μ = getϵμ( layer, matCol, wavenumber)
#     @assert isapprox(n^2, ϵ*μ, rtol=1e-3)
#
#     # Lecture 7B appears to take the conjugate of Kz
#     # Kz = conj( sqrt(I-Kx^2-Ky^2))
#     # This causes some tests to fail.  But seems to work when you also take the conjugate in the lambda calculation and
#     #TODO
#     return conj( ComplexF64[ sqrt(getk₀(kVectorSet)^2*conj(ϵ)*conj(μ) - kᵢ[X]^2 - kᵢ[Y]^2)  for kᵢ in kVectorSet.kᵢ] )  #TODO
#
#     # old:
#     # return -ComplexF64[ sqrt(getk₀(kVectorSet)^2*conj(ϵ)*conj(μ) - kᵢ[X]^2 - kᵢ[Y]^2)  for kᵢ in kVectorSet.kᵢ]
#
# end

function calc_ϵμ(layerDef::T1, matCol::MaterialCollection, kVectorSet::KVectorSet) where T1<:Union{UniformLayerDefinition, SemiInfiniteLayerDefinition}
    ϵ, μ = calc_ϵμ( getMaterial(matCol,layerDef.backgroundMaterialName), kVectorSet.wavenumber)
    return ϵ, μ
end

# Calling via layer definition.  For a uniform layer
function getϵμAtPosition( layerDef::T, positions::PositionGridXY{N}, materials::MaterialCollection, wavenumber::Wavenumber) where {N, T<:Union{SemiInfiniteLayerDefinition, UniformLayerDefinition}}

    mappedϵμ = map( position -> getϵμ( layerDef, materials, wavenumber), positions.positions )

    return mappedϵμ
end
# Calling via layer definition.  For a patterned layer
function getϵμAtPosition( layerDef::T, positionGrid::PositionGridXY{N}, materials::MaterialCollection, wavenumber::Wavenumber) where {N, T<:PatternedLayerDefinition}

    return map( position -> getϵμAtPosition( layerDef.layerPattern, position, materials, wavenumber), positionGrid.positions )
end



function getBoundaryLayer(layerStack, isTop)::SemiInfiniteLayerDefinition
    local layer
    if isTop
        layer = last(layerStack)
    else
        layer = first(layerStack)
    end
    return layer
end

# Get the index of the layerStack corresponding to the target Z position
# Why is this in lattice?
"""
    getLayerIndexFromZ(layerZpositions, zCoordinate)
Returns the index of the layer that contains the input zCoordinate, based on the Zpositions included.
"""
function getLayerIndexFromZ(layerZpositions, targetZ)
    # Layer N has z-layer positions N (bottom) and N+1 (top)
    if targetZ<first(layerZpositions)
        return 1
    elseif targetZ>last(layerZpositions)
        return length(layerZpositions)
    end
    return findfirst( layerZpositions .> targetZ ) - 1

end
