# Run all tests together.

include("UnitTests/TestSuite.jl")
include("IntegrationTests/IntegrationTest3.jl")
include("IntegrationTests/FresnelCoefficientTests.jl")
include("MaterialImport/MaterialDatabaseTests.jl")
