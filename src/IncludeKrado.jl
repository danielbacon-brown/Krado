using LinearAlgebra
using StaticArrays
using PyCall
using PyPlot
using Printf
using DelimitedFiles
using Interpolations
using YAML
# using TimerOutputs

# const global DIR = @__DIR__

export SimulationDefinition
export funcA
export getWavenumber


include("Miscellaneous\\Units.jl")
include("Miscellaneous\\Wavenumber.jl")
include("Miscellaneous\\MiscellaneousFunctions.jl")

# include("ModeField\\Mode.jl")



include("ScatteringMatrices\\QuadArray.jl")

include("Harmonics\\Lattice.jl")
include("Geometry\\PositionGrid.jl")
include("Harmonics\\Harmonics.jl")
include("Harmonics\\HarmonicsTruncation.jl")
include("Harmonics\\GvectorSet.jl")
include("Harmonics\\KVectorSet.jl")




include("Geometry\\Shapes.jl")
include("Materials\\Materials.jl")
include("Materials\\MaterialCollection.jl")
include("Materials\\MaterialDatabase.jl")
include("Geometry\\Solid.jl")

include("Layer\\LayerPattern.jl")
include("Layer\\GridLayerPattern.jl")
include("Layer\\LayerDefinition.jl")

include("Materials\\UserMaterials.jl")

include("ModeField\\ModeCoefficientSet.jl")
include("ModeField\\FieldSetStack.jl")
include("ModeField\\ElectricEigenvectors.jl")
include("ModeField\\FieldSetXYZ.jl")
include("ModeField\\FieldSetSP.jl")
include("ModeField\\FieldSetConversions.jl")

include("Layer\\LayerScatteringMatrix.jl")


include("Boundary\\SourceFields.jl")
include("Boundary\\BoundaryDefinition.jl")
include("Simulation\\Simulation.jl")

include("Analysis\\SimulationRunFunctions.jl")


include("Boundary\\BoundaryConditions.jl")

include("Simulation\\FreeSpaceParameters.jl")
include("Simulation\\DerivedParameters.jl")

include("ScatteringMatrices\\GlobalScatteringMatrix.jl")
include("ScatteringMatrices\\ScatteringMatrixCalculations.jl")

include("Analysis\\AnalysisFunctions.jl")
include("Analysis\\SimulationSweep.jl")

include("Plotting\\Plotting.jl")
include("Plotting\\CoordinateData.jl")
include("Plotting\\HarmonicsPlotting.jl")
include("Plotting\\LatticePlotting.jl")
include("Plotting\\LayerPlotting.jl")
include("Plotting\\Patch3D.jl")
include("Plotting\\LayerPatch3D.jl")
include("Plotting\\VectorPlotting.jl")
