# Integration test 3 with sweep
module IntegrationTest3
# Using RCWA-Benchmark-Data-3x3 by Raymond Rumpf from empossible.net

using Test



println()
println("Start:")

@testset "IntegrationTest3 Sweep" begin

println()
println("Beginning IntegrationTest3")

include("../../src/IncludeKrado.jl")


# Incident wavevector
λ₀ = 2*cm
wavenumber = WavenumberByλ₀(λ₀)

# Benchmark appears to use nonstandard rotation method.  Here, θ is azimuthal angle (rotation around z-axis) and ϕ is zenith angle (rotation around Y-axis).  ϕ rotation occurs first.
θ = 30 * degrees
ϕ = 60 * degrees
Eₚ = 0.70711
Eₛ = -0.70711im
inputAmplitudes = _2VectorComplex(Eₚ, Eₛ)
boundaryDefinition = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)

# Define material collection:
matCol = MaterialCollection()

# Reflection material:
bottomPermittivity = 2.0
bottomPermeability = 1.0
addMaterial!(matCol,"bottom material", Material( ConstantPermittivity(bottomPermittivity), ConstantPermeability(bottomPermeability) ) )

# Transmission material:
topPermittivity = 9.0
topPermeability = 1.0
addMaterial!(matCol,"top material", Material( ConstantPermittivity(topPermittivity), ConstantPermeability(topPermeability) ) )

# Device material:
devicePermittivity = 6.0
devicePermeability = 1.0
addMaterial!(matCol,"device material", Material( ConstantPermittivity(devicePermittivity), ConstantPermeability(devicePermeability) ) )

# Define lattice
Lx = 1.75 * cm
Ly = 1.5 * cm
U̅ = [Lx, 0]
V̅ = [0, Ly]
lattice = Lattice(U̅, V̅)

# Define layers.  Layer 1 contains triangle
layer1depth = 0.5 * cm  # thickness of layer 1
layer2depth = 0.3 * cm

# Define ayer 2:
layer2 = UniformLayerDefinition(layer2depth, "device material")

# Define layer 1:
# A grid of materialNames so that we can define a GridLayerPattern
w = 0.8*Ly  # side-length of triangle
Nx = 512
Ny = round(Int64,Nx*Ly/Lx) # Number of pixels in y-direction

dx = Lx/Nx
dy = Ly/Ny # Y-length of each pixel

materialNameGrid = String["device material" for ix=1:Nx, iy=1:Ny ]

h = 0.5*sqrt(3)*w # height of triangle
ny = round(Int64,h/dy) # Number of rows for triangle
ny1 = round(Int64,(Ny - ny)/2) # top y point of triangle
ny2 = ny1 + ny - 1 # bottom y edge of triangle
for ny = ny1:ny2
    f = (ny-ny1)/(ny2-ny1)
    nx = round(Int64,f*w/Lx*Nx)
    nx1 = 1 + floor(Int64,(Nx-nx)/2)
    nx2 = nx1 + nx
    materialNameGrid[nx1:nx2, ny] = String["bottom material" for ix=nx1:nx2]
end
layerPattern1 = GridLayerPattern( materialNameGrid )

# FFT real-space coordinates
Nx = 512
Ny = round(Int64, Nx*Ly/Lx)
layerDivisions = [Nx, Ny]
layer1 = PatternedLayerDefinition(layerDivisions, layer1depth, layerPattern1)

# Reflection and transmission layers
bottomLayer = SemiInfiniteLayerDefinition("bottom material")
topLayer = SemiInfiniteLayerDefinition("top material")

layerStack = [bottomLayer, layer1, layer2, topLayer]


# Define Harmonics.  Using 3x3-order only:
M,N = 1,1
harmonicsTruncation = HarmonicsTruncationByRectangle(M,N)

isForward = FORWARD
analysisDefinition = TransmittanceReflectanceAnalysisDefinition(isForward)




#### SWEEP OVER WAVELENGTH ####

# Create generation function for Simulation Definition.
# Parameters is an iterable
function doSimulationByWavelength(wavelength)

    local wavenumber = WavenumberByλ₀(wavelength)
    local boundaryDefinition = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)
    
    simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition)

    data = runSimulation(simulationDefinition)

    return data
end

# Set wavelengths in sweep
wavelengths = LinRange(1*cm, 3*cm, 3) 
sweepData = map(doSimulationByWavelength, wavelengths)

@test isapprox(sweepData[2].totalReflectance, 0.088768, rtol=1e-2)
@test isapprox(sweepData[2].totalTransmittance, 0.91123, rtol=1e-2)


#### SWEEP OVER GEOMETRY USING NAMEDTUPLE ####

# Create generation function for Simulation Definition.
# Parameters is a NamedTuple
function doSimulationByLatticeUsingNamedTuple(parameters)

    local U̅ = [parameters.Lx, 0]
    local V̅ = [0, parameters.Ly]
    local lattice = Lattice(U̅, V̅)

    simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition)

    data = runSimulation(simulationDefinition)

    return data
end

# Set lattice parameters in sweep
sweepParameters = (   (Lx = 1.5*cm, Ly = 1.25*cm),
                    (Lx = 1.75*cm, Ly = 1.5*cm),
                    (Lx = 2*cm, Ly = 1.75*cm), )      
sweepData = map(doSimulationByLatticeUsingNamedTuple, sweepParameters)

@test isapprox(sweepData[2].totalReflectance, 0.088768, rtol=1e-2)
@test isapprox(sweepData[2].totalTransmittance, 0.91123, rtol=1e-2)




#### SWEEP OVER GEOMETRY USING PARALLEL VECTORS ####

# Preferred method of sweeping over multiple parameters
function doSimulationByLattice(Lx, Ly)

    local U̅ = [Lx, 0]
    local V̅ = [0, Ly]
    local lattice = Lattice(U̅, V̅)

    simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition)

    data = runSimulation(simulationDefinition)

    return data
end

LxValues = LinRange(1.5*cm, 2*cm, 3)
LyValues = LinRange(1.25*cm, 1.75*cm, 3)
sweepData = map(doSimulationByLattice, LxValues, LyValues)

@test isapprox(sweepData[2].totalReflectance, 0.088768, rtol=1e-2)
@test isapprox(sweepData[2].totalTransmittance, 0.91123, rtol=1e-2)




#### 2D Sweep ###
# Run a wavelength sweep over wavelength and max harmonic orders M,N

function doSimulationByLattice(MN, λ₀)

    local harmonicsTruncation = HarmonicsTruncationByRectangle(MN,MN)

    local wavenumber = WavenumberByλ₀(λ₀)
    local boundaryDefinition = InputByOrderBoundaryDefinition(wavenumber, θ, ϕ, BOTTOM, inputAmplitudes)

    local simulationDefinition = SimulationDefinition(lattice, layerStack, harmonicsTruncation, boundaryDefinition, matCol, analysisDefinition)

    local data = runSimulation(simulationDefinition)

    return data
end

MNvalues = 0:2
λ₀Values = LinRange(1.5*cm, 2.5*cm, 3)
sweepData = map( prod -> doSimulationByLattice(prod[1], prod[2]), Iterators.product(MNvalues, λ₀Values))
sweepData = reshape(sweepData, (length(MNvalues),length(λ₀Values)) )
@test isapprox(sweepData[2,2].totalReflectance, 0.088768, rtol=1e-2)
@test isapprox(sweepData[2,2].totalTransmittance, 0.91123, rtol=1e-2)


end; # end of testset

end # end of module
