# A layer definition is a structure that describes how the convolution matrix should be calculated
# All have a thickness and a background material

export LayerDefinition, SemiInfiniteLayerDefinition, UniformLayerDefinition, PatternedLayerDefinition, getBoundaryLayer

abstract type LayerDefinition end

# # Used to determine how the points are positioned in the unit cell.
# abstract type GridAlignment end
# mutable struct LeftAlignment <: GridAlignment end
# mutable struct CenterAlignment <: GridAlignment end
# const global LEFTALIGNMENT = LeftAlignment()
# const global CENTERALIGNMENT = CenterAlignment()


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
    # gridAlignment::GridAlignment

    function PatternedLayerDefinition(numDivisions::_2VectorInt, thickness::Real, layerPattern::AbstractLayerPattern)
        return new(convert(_2VectorInt,numDivisions), convert(Float64,thickness), layerPattern)
    end
    # function PatternedLayerDefinition(numDivisions::_2VectorInt, thickness::Real, layerPattern::AbstractLayerPattern, gridAlignment::GridAlignment)
    #     return new(convert(_2VectorInt,numDivisions), convert(Float64,thickness), layerPattern, gridAlignment)
    # end
end
# function PatternedLayerDefinition(numDivisions::_2VectorInt, thickness::Real, layerPattern::AbstractLayerPattern; gridAlignment::GridAlignment = CENTERALIGNMENT )
#     return new(convert(_2VectorInt,numDivisions), convert(Float64,thickness), layerPattern, gridAlignment)
# end
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








function getPositionϵμGrids( layerDef::PatternedLayerDefinition, lattice::Lattice, matCol::MaterialCollection, wavenumber::Wavenumber)
    # positionGrid = PositionGridXYleftAligned(lattice, layerDef.numDivisions)
    positionGrid = PositionGridXY(lattice, layerDef.numDivisions)
    ϵGrid = getϵAtPosition( layerDef.layerPattern, positionGrid, matCol, wavenumber)
    μGrid = getμAtPosition( layerDef.layerPattern, positionGrid, matCol, wavenumber)
    return positionGrid, ϵGrid, μGrid
end








# Returns the position grid Array{_2VectorFloat(X,Y),2}
# function calcUniformGridPositions(lattice::Lattice, layerDef::PatternedLayerDefinition)
#     return calcUniformGridPositions(lattice, layerDef.numDivisions)
# end


#shorthand for gettign the refractive index of a uniform layer
function getn(layer::layerT, matCol::MaterialCollection, wavenumber::Wavenumber) where layerT<:Union{SemiInfiniteLayerDefinition, UniformLayerDefinition}
    return convert_ϵ2n(calc_ϵ( getMaterial(matCol,layer.backgroundMaterialName), wavenumber))
end

function getϵμ(layer::layerT, matCol::MaterialCollection, wavenumber::Wavenumber) where layerT<:Union{SemiInfiniteLayerDefinition, UniformLayerDefinition}
    return calc_ϵμ( getMaterial(matCol,layer.backgroundMaterialName), wavenumber)
end

function calckzBottom(kVectorSet::KVectorSet, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection, wavenumber::Wavenumber)
    ϵ,μ = getϵμ( layer, matCol, wavenumber)
    return ComplexF64[ -sqrt.( conj(ϵ)*conj(μ) - kᵢ[X]^2 - kᵢ[Y]^2)  for kᵢ in kVectorSet.kᵢNorm]  # REMOVING CONJ from Kzᵦ.  Either the lecture or the benchmark is wrong.
end

function calckzTop(kVectorSet::KVectorSet, layer::SemiInfiniteLayerDefinition, matCol::MaterialCollection, wavenumber::Wavenumber)
    ϵ,μ = getϵμ( layer, matCol, wavenumber)
    return ComplexF64[ sqrt.( conj(ϵ)*conj(μ) - kᵢ[X]^2 - kᵢ[Y]^2)  for kᵢ in kVectorSet.kᵢNorm]  # REMOVING CONJ from Kzᵦ.  Either the lecture or the benchmark is wrong.
end

# function calc_ϵμ(layerDef::T1, matCol::MaterialCollection, kVectorSet::KVectorSet) where T1<:Union{UniformLayerDefinition, SemiInfiniteLayerDefinition}
function calc_ϵμ(layerDef::T1, matCol::MaterialCollection, wavenumber::Wavenumber) where T1<:Union{UniformLayerDefinition, SemiInfiniteLayerDefinition}
    ϵ, μ = calc_ϵμ( getMaterial(matCol,layerDef.backgroundMaterialName), wavenumber)
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
