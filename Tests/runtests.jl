# Run all tests together.

include("UnitTests/TestSuite.jl")
include("IntegrationTests/FresnelCoefficientTests.jl")
include("IntegrationTests/IntegrationTest3.jl")
include("IntegrationTests/FieldSetConversionTests.jl")
include("MaterialImport/MaterialDatabaseTests.jl")
