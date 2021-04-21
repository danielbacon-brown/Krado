using LinearAlgebra
using StaticArrays
using PyCall
using PyPlot
# using Plots
using Printf
using DelimitedFiles
using Interpolations
using YAML
# using TimerOutputs
# using pyimport

# const global DIR = @__DIR__

export SimulationDefinition
export funcA
export getWavenumber


include("Miscellaneous\\Units.jl")
include("Miscellaneous\\Wavenumber.jl")
include("Miscellaneous\\MiscellaneousFunctions.jl")

include("ModeField\\Mode.jl")



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

include("Materials\\FavoriteMaterials.jl")
# include("..\\UserMaterials\\UserMaterialImporters.jl")  # This has been moved outside the package.  The user must now include the functions according to the location of their files.

include("ModeField\\ModeCoefficientSet.jl")
include("ModeField\\ModeFieldSet.jl")
include("ModeField\\ElectricEigenvectors.jl")
include("ModeField\\ModeField3Set.jl")
include("ModeField\\ModeFieldSPSet.jl")

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