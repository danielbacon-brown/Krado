# Run all tests together.

include("UnitTests/TestSuite.jl")
include("IntegrationTests/IntegrationTest3.jl")
# include("TestsOnExternals/AngledReflectionFromSi.jl")
include("IntegrationTests/SandPtransmission.jl")
# include("TestsOnExternals/MaterialDatabaseTests.jl")
