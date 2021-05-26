# Run all tests together.
module RunTests

include("UnitTests/TestSuite.jl")
include("IntegrationTests/FresnelCoefficientTests.jl")
include("IntegrationTests/IntegrationTest3.jl")
include("IntegrationTests/FieldSetConversionTests.jl")
include("MaterialImport/MaterialDatabaseTests.jl")

end;
