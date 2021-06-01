module TestModule

using Test

println()

# relDir = "Tests/UnitTests/"
relDir = ""

# include("../../src/Krado.jl")
# using .Krado
include("../../src/IncludeKrado.jl")

@testset "Unit Tests" begin

include(relDir * "LayerDefinitionTests.jl")

# include(relDir * "MaterialDatabaseTests.jl")  # Not including this in test suite because it relies on external files.

# include(relDir * "ModeTests.jl")  # NO LONGER USING "MODE"

include(relDir * "UnitsTests.jl")

include(relDir * "ShapesTests.jl")

include(relDir * "LayerSetTests.jl")

include(relDir * "ScatteringMatrixTests.jl")

include(relDir * "BoundaryConditionsTests.jl")

include(relDir * "KVectorSetTests.jl")

include(relDir * "MiscellaneousTests.jl")

include(relDir * "NestedArraysTests.jl")

include(relDir * "LatticeTests.jl")

include(relDir * "HarmonicsSetTests.jl")
include(relDir * "RectangularHarmonicsTests.jl")
include(relDir * "PowerHarmonicsTests.jl")
include(relDir * "GvectorsTests.jl")

include(relDir * "ShapesTests.jl")

include(relDir * "MaterialsTests.jl")


include(relDir * "LayerPatternTests.jl")


end;

end
