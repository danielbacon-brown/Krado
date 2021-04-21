module IntegrationTest2
# Using RCWA-Benchmark-Data-3x3 by Raymond Rumpf from empossible.net

using Test
using LinearAlgebra
using StaticArrays
using Printf


println()
println()

include("IntegrationTest2BenchmarkData.jl")

@testset "IntegrationTest2" begin

    include("../IncludeKrado.jl")

    # Incident wavevector
    λ₀ = 2*cm
    θ = 0 * degrees
    phi = 0 * degrees
    Tₚpower = 1
    Tₘpower = 0

    # Define materials:
    matCol = MaterialCollection()

    # Reflection material:
    refl_permittivity = 2.0
    refl_permeability = 1.0
    addMaterial!(matCol,"refl_material", Material( ConstantPermittivity(refl_permittivity), ConstantPermeability(refl_permeability) ) )

    # Transmission material:
    tran_permittivity = 9.0
    tran_permeability = 1.0
    addMaterial!(matCol,"tran_material", Material( ConstantPermittivity(tran_permittivity), ConstantPermeability(tran_permeability) ) )

    # Device material:
    device_permittivity = 6.0
    device_permeability = 1.0
    addMaterial!(matCol,"device_material", Material( ConstantPermittivity(device_permittivity), ConstantPermeability(device_permeability) ) )



    # Pitch
    Lx = 1.75 * cm
    Ly = 1.5 * cm
    U̅ = [Lx, 0]
    V̅ = [0, Ly]
    lattice = Lattice(U̅, V̅)

    # Harmonics definition
    # Using 3x3-order only:
    M = 1
    N = 1
    harmonicsDef = HarmonicsTruncationByRectangle(M,N)
    harmonicsSet = calcHarmonicsSet(harmonicsDef)
    @test harmonicsSet.mnᵢ == [_2VectorInt(-1, -1), _2VectorInt(0, -1), _2VectorInt(1, -1), _2VectorInt(-1,0), _2VectorInt(0, 0), _2VectorInt(1,0), _2VectorInt(-1, 1), _2VectorInt(0,1), _2VectorInt(1, 1),]
    @test harmonicsSet.indᵢ_mn == Dict(_2VectorInt(-1, -1)=>1, _2VectorInt(0,-1)=>2, _2VectorInt(1, -1)=>3, _2VectorInt(-1,0)=>4, _2VectorInt(0, 0)=>5, _2VectorInt(1,0)=>6, _2VectorInt(-1, 1)=>7, _2VectorInt(0,1)=>8, _2VectorInt(1, 1)=>9)

    #K-vectors and G-vectors
    Gvectors = GvectorSet(harmonicsSet,lattice)


    # Define layers.  Layer 1 contains triangle
    layer1_d = 0.5 * cm  # thickness of layer 1
    layer2_d = 0.3 * cm

    # Define unpatterned layer:
    layer2 = UniformLayerDefinition(layer2_d, "device_material")


    # Define shape: equilateral triangle pointed toward +y.

    # Reference style:
    # Define a grid of materialNames so that we can define a GridLayerPattern
    w = 0.8*Ly  # width of triangle
    Nx = 512
    Ny = round(Int64,Nx*Ly/Lx) # Number of pixels in y-direction

    dx = Lx/Nx
    dy = Ly/Ny # Y-length of each pixel

    materialNameGrid = String["device_material" for ix=1:Nx, iy=1:Ny ]

    h = 0.5*sqrt(3)*w # height of triangle
    ny = round(Int64,h/dy) # Number of rows for triangle
    ny1 = round(Int64,(Ny - ny)/2) # top y point of triangle
    ny2 = ny1 + ny - 1 # bottom y edge of triangle
    for ny = ny1:ny2
        f = (ny-ny1)/(ny2-ny1)
        nx = round(Int64,f*w/Lx*Nx)
        nx1 = 1 + floor(Int64,(Nx-nx)/2)
        nx2 = nx1 + nx
        materialNameGrid[nx1:nx2, ny] = String["refl_material" for ix=nx1:nx2]
    end
    layerPattern1 = GridLayerPattern( materialNameGrid )

    # simple style:
    # triangle_length = 0.8 * Ly
    # point1 = [triangle_length/2, -(triangle_length/2) / tand(60)]
    # point2 = [-triangle_length/2, -(triangle_length/2) / tand(60)]
    # point3 = [0, (triangle_length/2) / sind(60)]
    # triangle_shape = Polygon([Lx/2,Ly/2],[point1, point2, point3])
    # triangle_solid = Solid( triangle_shape, "refl_material" )
    # # Define solids calculator
    # layerPattern1 = LayerPattern( triangle_solid, "device_material")


    # FFT real-space coordinates
    # Nx = 20
    Nx = 512
    Ny = round(Int64, Nx*Ly/Lx)
    # @test Ny == 439
    layerDivisions = [Nx, Ny]
    layer1 = PatternedLayerDefinition(layerDivisions, layer1_d, layerPattern1)

    # Reflection and transmission layers

    layerReflection = SemiInfiniteLayerDefinition("refl_material")
    layerTransmission = SemiInfiniteLayerDefinition("tran_material")


    # PLOT positions and patterned grid:
    # posGrid = calcUniformGridPositions(lattice, layer1)
    # plotPositionScatter(layer1, lattice)
    # ϵGrid, μGrid = calcϵμGrids( lattice, layer1, matDict, λ₀)
    # plotGrid(ϵGrid, "ϵ' of layer 1", "ϵ'' of layer 1", "X", "Y")
    # plotGridByScatter(real(ϵGrid), posGrid, "ϵ of layer 1", "X", "Y")


    # Calculate convolution matrices:
    Cϵᵢⱼ1, Cμᵢⱼ1 = calcConvolutionMatrices( layer1, lattice, Gvectors, matCol, λ₀ )
    # @test Cμᵢⱼ1[1,1] ≈ Complex(1,0)
    @test Cμᵢⱼ1 ≈ Array{ComplexF64,2}(I,(9,9))

    #Original orientation:
    Cϵᵢⱼ1_benchmark = ComplexF64[
        5.0449+0.0im	0.636-0.002im	-0.1402+0.0008im	0.3671+0.4094im	-0.3056-0.1696im	0.1402-0.1531im	0.2044-0.0362im	-0.0686+0.0227im	-0.0727+0.0443im;
        0.636+0.002im	5.0449+0.0im	0.636-0.002im	-0.3045-0.1715im	0.3671+0.4094im	-0.3056-0.1696im	-0.0687+0.0222im	0.2044-0.0362im	-0.0686+0.0227im;
        -0.1402-0.0008im	0.636+0.002im	5.0449+0.0im	0.1421-0.1514im	-0.3045-0.1715im	0.3671+0.4094im	-0.0733+0.0433im	-0.0687+0.0222im	0.2044-0.0362im;
        0.3671-0.4094im	-0.3045+0.1715im	0.1421+0.1514im	5.0449+0.0im	0.636-0.002im	-0.1402+0.0008im	0.3671+0.4094im	-0.3056-0.1696im	0.1402-0.1531im;
        -0.3056+0.1696im	0.3671-0.4094im	-0.3045+0.1715im	0.636+0.002im	5.0449+0.0im	0.636-0.002im	-0.3045-0.1715im	0.3671+0.4094im	-0.3056-0.1696im;
        0.1402+0.1531im	-0.3056+0.1696im	0.3671-0.4094im	-0.1402-0.0008im	0.636+0.002im	5.0449+0.0im	0.1421-0.1514im	-0.3045-0.1715im	0.3671+0.4094im;
        0.2044+0.0362im	-0.0687-0.0222im	-0.0733-0.0433im	0.3671-0.4094im	-0.3045+0.1715im	0.1421+0.1514im	5.0449+0.0im	0.636-0.002im	-0.1402+0.0008im;
        -0.0686-0.0227im	0.2044+0.0362im	-0.0687-0.0222im	-0.3056+0.1696im	0.3671-0.4094im	-0.3045+0.1715im	0.636+0.002im	5.0449+0.0im	0.636-0.002im;
        -0.0727-0.0443im	-0.0686-0.0227im	0.2044+0.0362im	0.1402+0.1531im	-0.3056+0.1696im	0.3671-0.4094im	-0.1402-0.0008im	0.636+0.002im	5.0449+0.0im]

    @test isapprox(Cϵᵢⱼ1, Cϵᵢⱼ1_benchmark, rtol=1e-1)


    Cϵᵢⱼ2, Cμᵢⱼ2 = calcConvolutionMatrices( layer2, lattice, Gvectors, matCol, λ₀ )
    @test Cϵᵢⱼ2[1,1] ≈ Complex(6,0)
    @test Cμᵢⱼ2[1,1] ≈ Complex(1,0)
    @test Cϵᵢⱼ2 ≈ 6*Array{ComplexF64,2}(I,(9,9))
    @test Cμᵢⱼ2 ≈ Array{ComplexF64,2}(I,(9,9))

    #Calculate inverse convolution matrices for each layer
    Cϵᵢⱼ⁻¹1 = inv(Cϵᵢⱼ1)
    Cμᵢⱼ⁻¹1 = inv(Cμᵢⱼ1)
    Cϵᵢⱼ⁻¹2 = inv(Cϵᵢⱼ2)
    Cμᵢⱼ⁻¹2 = inv(Cμᵢⱼ2)





    # Wave vector expansion:
    nᵣ = convert_ϵ2n(calc_ϵ( getMaterial(matCol,"refl_material"), λ₀))
    nₜ = convert_ϵ2n(calc_ϵ( getMaterial(matCol,"tran_material"), λ₀))
    @test isapprox( nᵣ, 1.4142, rtol = 1e-3)
    @test isapprox( nₜ, 3, rtol = 1e-3)
    @test isapprox( λ₀2k₀(λ₀), 314.1593, rtol = 1e-3)

    # Create the incident wave-vector
    k₀ = λ₀2k₀(λ₀)
    k = normalIncidenceKvector(k₀)
    @test norm(k) ≈ k₀
    kVectorSet = createKVectorSet(k, Gvectors)
    
    # Calculate the z-components of k-vector in the transmitted and reflected layer
    kzᵣ = Diagonal( ComplexF64[ -1*conj(sqrt((kVectorSet.k₀*nᵣ)^2 - kᵢ[X]^2 - kᵢ[Y]^2))  for kᵢ in kVectorSet.kᵢ] ) 
    @test isapprox(kzᵣ/kVectorSet.k₀, kzᵣbenchmark, rtol=1e-3)
    kzₜ = Diagonal( ComplexF64[ conj(sqrt((kVectorSet.k₀*nₜ)^2 - kᵢ[X]^2 - kᵢ[Y]^2))  for kᵢ in kVectorSet.kᵢ] ) 
    @test isapprox(kzₜ/kVectorSet.k₀, kzₜbenchmark, rtol=1e-3)
    
    kzᵣ = Diagonal( calckzᵣ(kVectorSet, layerReflection, matCol, λ₀) ) 
    @test isapprox(kzᵣ/kVectorSet.k₀, kzᵣbenchmark, rtol=1e-3)
    kzₜ = Diagonal( calckzₜ(kVectorSet, layerTransmission, matCol, λ₀) ) 
    @test isapprox(kzₜ/kVectorSet.k₀, kzₜbenchmark, rtol=1e-3)
    
    
    
    # print(kVectorSet.Kx)

    # DEFINE INPUT MODE

    # Using TE as parallel to y
    âTE = _2VectorComplex(0,1) # Polarization amplitude of TE vector
    âTM = _2VectorComplex(1,0) # Polarization amplitude of TM vector
    inputModeTE = Mode(k,âTE)
    inputModeTM = Mode(k,âTM)

    # DEFINE INPUT BOUNDARY CONDITIONS
    inputOrder = _2VectorInt(0,0)
    boundaryConditionsTE = InputByOrderBoundaryConditions(inputModeTE)
    boundaryConditionsTM = InputByOrderBoundaryConditions(inputModeTM, inputOrder)

    @test k[X] == 0
    @test k[Y] == 0

    #STEP 5: Calc eigenmodes for free space
    # same for all layers
    KzNorm = calcKzForUnpatternedEigenmode(kVectorSet)
    @test isapprox(KzNorm, KzNormBenchmark, rtol=1e-3)

    Q = calcQForUnpatternedEigenmode( kVectorSet )
    @test isapprox(Q, Q_benchmark, rtol=1e-3)

    W₀ = calcW₀( numHarmonics(kVectorSet) )
    @test W₀ ≈ Array{ComplexF64,2}(I,(9*2,9*2))

    λeigenvalues = calcλForUnpatternedEigenMode(KzNorm)
    @test isapprox(λeigenvalues, λeigenvalues_benchmark, rtol=1e-3)

    V₀ = calcEigenmodesFromQλ(Q,λeigenvalues)
    @test isapprox(V₀, V₀benchmark, rtol=1e-3)


    # V₀ = calcEigenmodesForUnpatterned(kVectorSet)
    V₀ = calcV₀(kVectorSet)
    @test isapprox(V₀, V₀benchmark, rtol=1e-3)



    #STEP 6: Initialize global S
    Sglobal = initializeGlobalS(numHarmonics(harmonicsSet))
    SglobalBenchmark11 = zeros(ComplexF64,(9*2,9*2))
    SglobalBenchmark12 = Array{ComplexF64,2}(I,(9*2,9*2))
    SglobalBenchmark21 = Array{ComplexF64,2}(I,(9*2,9*2))
    SglobalBenchmark22 = zeros(ComplexF64,(9*2,9*2))
    SglobalBenchmark = vcat( hcat(SglobalBenchmark11, SglobalBenchmark12),
                            hcat(SglobalBenchmark21, SglobalBenchmark22) )
    @test Sglobal ≈ SglobalBenchmark

    # STEP 7: MAIN LOOP

    # Layer 1: Patterned
    P₁ = calcPmatrixPatterned(kVectorSet, Cϵᵢⱼ1, Cϵᵢⱼ⁻¹1, Cμᵢⱼ1, Cμᵢⱼ⁻¹1)
    @test isapprox(P₁, P₁benchmark, rtol=1e-3)

    Q₁ = calcQmatrixPatterned(kVectorSet, Cϵᵢⱼ1, Cϵᵢⱼ⁻¹1, Cμᵢⱼ1, Cμᵢⱼ⁻¹1)
    @test isapprox(Q₁, Q₁benchmark, rtol=1e-3)

    Ω²₁ = calcΩ²(P₁,Q₁)
    @test isapprox(Ω²₁, Ω²₁benchmark, rtol=1e-3)

                        
    # NOTE: The eigenvalue decomposition is not unique.  Should not expect these tests to pass
    W₁, λ₁ = calcWᵢλᵢ(Ω²₁)
    @test Ω²₁ * W₁ ≈ W₁ * (λ₁.^2)  # Confirming that it is an eigendecomposition
    @test isapprox(Ω²₁benchmark * W₁benchmark, W₁benchmark * (λ₁benchmark.^2), rtol=1e-3)
    @test isapprox(Ω²₁benchmark * W₁benchmark, W₁benchmark * (λ₁benchmark^2), rtol=1e-3)



    V₁ = calcMagneticEigenvectorsFromQWλ(Q₁, W₁, λ₁)
    # @test isapprox(V₁, V₁benchmark, rtol=1e-3)

    A₁, B₁ = calcAB(W₁,W₀,V₁,V₀)
    X₁ = calcX(λ₁, kVectorSet.k₀, layer1.thickness)

    _1, _2 = getQuadrantSlices(numHarmonics(kVectorSet))

    S₁ = calcScatteringMatrix_ABX(A₁,B₁,X₁)

    @test isapprox(S₁[_1,_1], S1₁₁benchmark, rtol=1e-3)
    @test isapprox(S₁[_1,_2], S1₁₂benchmark, rtol=1e-3)
    @test isapprox(S₁[_2,_1], S1₂₁benchmark, rtol=1e-3)
    @test isapprox(S₁[_2,_2], S1₂₂benchmark, rtol=1e-3)

    # sugary:
    S₁ = calcScatteringMatrix(layer1, matCol, kVectorSet, Gvectors, lattice )
    @test isapprox(S₁[_1,_1], S1₁₁benchmark, rtol=1e-3)
    @test isapprox(S₁[_1,_2], S1₁₂benchmark, rtol=1e-3)
    @test isapprox(S₁[_2,_1], S1₂₁benchmark, rtol=1e-3)
    @test isapprox(S₁[_2,_2], S1₂₂benchmark, rtol=1e-3)


    # Layer 2: Unpatterned
    mat₂str = layer2.backgroundMaterialName
    mat₂ = getMaterial(matCol,mat₂str)
    ϵ₂, μ₂ = calc_ϵμ(mat₂,λ₀)
    P₂ = calcPmatrixUnpatterned(kVectorSet, ϵ₂, μ₂ )
    @test isapprox(P₂, P₂benchmark, rtol=1e-3) 

    Q₂ = calcQmatrixUnpatterned(P₂, ϵ₂, μ₂)
    @test isapprox(Q₂, Q₂benchmark, rtol=1e-3) 


    # Sugary version.  Method specific to uniform layer
    P₂, Q₂ = calcPQmatrix(layer2, kVectorSet, matCol)
    @test isapprox(P₂, P₂benchmark, rtol=1e-3) 
    @test isapprox(Q₂, Q₂benchmark, rtol=1e-3) 

    # calcΩ² same for each layer
    Ω²₂ = calcΩ²(P₂,Q₂)
    @test isapprox(Ω²₂, Ω²₂benchmark, rtol=1e-3) 

    W₂, λ₂ = calcWᵢλᵢ(Ω²₂)
    @test isapprox(Ω²₂benchmark * W₂benchmark, W₂benchmark * (λ₂benchmark.^2), rtol=1e-3)
    # @test isapprox(W₂, W₂benchmark, rtol=1e-3) 
    # @test isapprox(λ₂, λ₂benchmark, rtol=1e-3) 

    V₂ = calcMagneticEigenvectorsFromQWλ(Q₂,W₂,λ₂)
    # @test isapprox(V₂, V₂benchmark, rtol=1e-3)

    # sugary:
    V₂ = calcEigenmodesForUniformLayer(kVectorSet, layer2, matCol)
    # @test isapprox(V₂, V₂benchmark, rtol=1e-3)


    # Common components of scattering matrix
    A₂ = calcA(W₂, W₀, V₂, V₀)
    # @test isapprox(A₂, A₂benchmark, rtol=1e-3)

    B₂ = calcB(W₂, W₀, V₂, V₀)
    # @test isapprox(B₂, B₂benchmark, rtol=1e-3)

    X₂ = calcX(λ₂, kVectorSet.k₀, layer2.thickness)
    # @test isapprox(X₂, X₂benchmark, rtol=1e-3)


    S₂ = calcScatteringMatrix_ABX(A₂, B₂, X₂)
    @test isapprox(S₂[_1,_1], S2₁₁benchmark, rtol=1e-3)
    @test isapprox(S₂[_1,_2], S2₁₂benchmark, rtol=1e-3)
    @test isapprox(S₂[_2,_1], S2₂₁benchmark, rtol=1e-3)
    @test isapprox(S₂[_2,_2], S2₂₂benchmark, rtol=1e-3)

    
    S₂ = calcScatteringMatrix(layer2, matCol, kVectorSet)
    @test isapprox(S₂[_1,_1], S2₁₁benchmark, rtol=1e-3)
    @test isapprox(S₂[_1,_2], S2₁₂benchmark, rtol=1e-3)
    @test isapprox(S₂[_2,_1], S2₂₁benchmark, rtol=1e-3)
    @test isapprox(S₂[_2,_2], S2₂₂benchmark, rtol=1e-3)



    # STEP 8: Reflection side scattering matrix
    Pᵣ, Qᵣ = calcPQmatrix(layerReflection, kVectorSet, matCol)
    @test isapprox(Qᵣ, Qᵣbenchmark, rtol=1e-3) 
    Ω²ᵣ = calcΩ²(Pᵣ,Qᵣ)
    # @test isapprox(Ω²ᵣ, Ω²ᵣbenchmark, rtol=1e-3)
    
    λᵣ = Array(vcat( hcat(-1im*kzᵣ, zeros(ComplexF64,size(kzᵣ)) ),
               hcat(zeros(ComplexF64,size(kzᵣ )), -1im*kzᵣ) ) / kVectorSet.k₀)
    @test isapprox(λᵣ, λᵣbenchmark, rtol=1e-3)
    λᵣ = calcλreflection(kzᵣ, kVectorSet.k₀)
    @test isapprox(λᵣ, λᵣbenchmark, rtol=1e-3)
    Wᵣ = W₀
        

    # Vᵣ = calcEigenmodesFromQλ(Qᵣ,λᵣ)
    Vᵣ = calcMagneticEigenvectorsFromQWλ(Qᵣ,Wᵣ,λᵣ)
    @test isapprox(Vᵣ,Vᵣbenchmark,rtol=1e-3)

    Aᵣ = calcA_SemiInfinite(Wᵣ, W₀, Vᵣ, V₀)
    @test isapprox(Aᵣ, Aᵣbenchmark, rtol=1e-3)

    Bᵣ = calcB_SemiInfinite(Wᵣ, W₀, Vᵣ, V₀)
    @test isapprox(Bᵣ, Bᵣbenchmark, rtol=1e-3)

    # #sugary
    Aᵣ, Bᵣ = calcAB_SemiInfinite(Wᵣ, W₀, Vᵣ, V₀)
    @test isapprox(Aᵣ, Aᵣbenchmark, rtol=1e-3)
    @test isapprox(Bᵣ, Bᵣbenchmark, rtol=1e-3)

    Sᵣ = calcScatteringMatrixReflection_AB(Aᵣ,Bᵣ)
    @test isapprox(Sᵣ[_1,_1], SR₁₁benchmark, rtol=1e-3)
    @test isapprox(Sᵣ[_1,_2], SR₁₂benchmark, rtol=1e-3)
    @test isapprox(Sᵣ[_2,_1], SR₂₁benchmark, rtol=1e-3)
    @test isapprox(Sᵣ[_2,_2], SR₂₂benchmark, rtol=1e-3)


    Sᵣ = calcScatteringMatrixReflection(layerReflection, matCol, kVectorSet)
    @test isapprox(Sᵣ[_1,_1], SR₁₁benchmark, rtol=1e-3)
    @test isapprox(Sᵣ[_1,_2], SR₁₂benchmark, rtol=1e-3)
    @test isapprox(Sᵣ[_2,_1], SR₂₁benchmark, rtol=1e-3)
    @test isapprox(Sᵣ[_2,_2], SR₂₂benchmark, rtol=1e-3)


    # STEP 9: Transmission side scattering matrix
    Pₜ, Qₜ = calcPQmatrix(layerTransmission, kVectorSet, matCol)
    Ω²ₜ = calcΩ²(Pₜ,Qₜ)

    λₜ = calcλtransmission(kzₜ, kVectorSet.k₀)
    @test isapprox(λₜ,λₜbenchmark,rtol=1e-3)
    Wₜ = W₀

    Vₜ = calcMagneticEigenvectorsFromQWλ(Qₜ,Wₜ,λₜ)

    Aₜ = calcA_SemiInfinite(Wₜ, W₀, Vₜ, V₀)
    Bₜ = calcB_SemiInfinite(Wₜ, W₀, Vₜ, V₀)

    Sₜ = calcScatteringMatrixTransmission_AB(Aₜ,Bₜ)
    @test isapprox(Sₜ[_1,_1], ST₁₁benchmark, rtol=1e-3)
    @test isapprox(Sₜ[_1,_2], ST₁₂benchmark, rtol=1e-3)
    @test isapprox(Sₜ[_2,_1], ST₂₁benchmark, rtol=1e-3)
    @test isapprox(Sₜ[_2,_2], ST₂₂benchmark, rtol=1e-3)


    Sₜ = calcScatteringMatrixTransmission(layerTransmission, matCol, kVectorSet)


    # STEP 10: Global scattering matrix:
    Sdevice = S₁⊗S₂
    Sglobal = Sᵣ⊗Sdevice⊗Sₜ
    @test isapprox(Sglobal[_1,_1], SG₁₁benchmark,rtol=1e-3)
    @test isapprox(Sglobal[_1,_2], SG₁₂benchmark,rtol=1e-3)
    @test isapprox(Sglobal[_2,_1], SG₂₁benchmark,rtol=1e-3)
    @test isapprox(Sglobal[_2,_2], SG₂₂benchmark,rtol=1e-3)
    

    # Put it in terms of a device stack:
    layerStack = [layerReflection, layer1, layer2, layerTransmission]
    Sglobal = calcScatteringMatrix(layerStack, matCol, kVectorSet, Gvectors, lattice)
    @test isapprox(Sglobal[_1,_1], SG₁₁benchmark,rtol=1e-3)
    @test isapprox(Sglobal[_1,_2], SG₁₂benchmark,rtol=1e-3)
    @test isapprox(Sglobal[_2,_1], SG₂₁benchmark,rtol=1e-3)
    @test isapprox(Sglobal[_2,_2], SG₂₂benchmark,rtol=1e-3)


    # STEP 11: FIELDS
    ŝ, p̂ = calcŝp̂(getRealkXYZ(inputModeTE))
    @test ŝ == [1,0,0]
    @test p̂ == [0,1,0]
    ŝ, p̂ = calcŝp̂(getRealkXYZ(inputModeTM))
    @test ŝ == [1,0,0]
    @test p̂ == [0,1,0]

    # NOTE: The benchmark uses TE as parallel to p̂
    𝐏TE = fieldSPtoFieldXYZ(inputModeTE)
    𝐏TM = fieldSPtoFieldXYZ(inputModeTM)
    @test 𝐏TE ≈ [0,1,0]
    @test 𝐏TM ≈ [1,0,0]


    # Source fields:
    sourceFields1TE, sourceFields2TE = calcSourceFields(boundaryConditionsTE, harmonicsSet, kVectorSet)

    
    sourceFields1TEbenchmark = [0;0;0;0;0;0;0;0;0;
                                    0;0;0;0;1;0;0;0;0]
    sourceFields2TEbenchmark = [0;0;0;0;0;0;0;0;0;
                                    0;0;0;0;0;0;0;0;0]

    @test isapprox(sourceFields1TE, sourceFields1TEbenchmark,rtol=1e-3 )
    @test isapprox(sourceFields2TE, sourceFields2TEbenchmark,rtol=1e-3 )


    sourceFields1TM, sourceFields2TM = calcSourceFields(boundaryConditionsTM, harmonicsSet, kVectorSet)


    sourceModeCoeff1TE = modeField2Coeff(sourceFields1TE, Wᵣ)
    
    sourceModeCoeff1TEbenchmark = [0;0;0;0;0;0;0;0;0;
                                    0;0;0;0;1;0;0;0;0]

    sourceModeCoeff2TE = modeField2Coeff(sourceFields2TE, Wₜ)
    @test isapprox( sourceModeCoeff1TE, sourceModeCoeff1TEbenchmark, rtol=1e-3)
    
    S₁₁ = Sglobal[_1,_1]
    S₁₂ = Sglobal[_1,_2]
    S₂₁ = Sglobal[_2,_1]
    S₂₂ = Sglobal[_2,_2]

    # Compute transmission and reflection mode coefficients ("cref" in benchmark)
    reflModeCoeffTE = S₁₁*sourceModeCoeff1TE
    transModeCoeffTE = S₂₁*sourceModeCoeff1TE
    reflModeCoeffTEbenchmark = [
        -0.0094 - 0.0369im;	0.0000 + 0.0000im;	0.0092 + 0.0369im;	0.0012 - 0.0042im;	0.0000 + 0.0000im;	-0.0012 + 0.0042im;	0.0374 + 0.0158im;0.0000 - 0.0000im;	-0.0373 - 0.0160im;	0.0432 + 0.0301im;	-0.0265 - 0.0191im;	0.0430 + 0.0304im;-0.1028 + 0.0241im;-0.2706 - 0.0715im;-0.1029 + 0.0235im;	0.0378 - 0.0313im;	-0.0199 + 0.0283im;	0.0380 - 0.0310im]
    transModeCoeffTEbenchmark = [
            0.0268 + 0.0588im;0.0000 + 0.0000im;-0.0264 - 0.0589im;-0.0134 + 0.0036im;0.0000 + 0.0000im;0.0134 - 0.0035im;-0.0569 - 0.0036im;-0.0000 + 0.0000im;0.0569 + 0.0039im;-0.0263 - 0.0824im;0.0699 + 0.1017im;-0.0258 - 0.0825im;0.1382 - 0.0589im;0.5145 + 0.2249im;0.1386 - 0.0581im;-0.0896 - 0.0165im;0.0809 - 0.0620im;-0.0895 - 0.0171im ]
    @test isapprox(reflModeCoeffTE, reflModeCoeffTEbenchmark, rtol=1e-2)
    @test isapprox(transModeCoeffTE, transModeCoeffTEbenchmark, rtol=1e-2)



    # simpler:
    reflModeCoeffTE, transModeCoeffTE = propagateModeCoeff(Sglobal, sourceModeCoeff1TE, sourceModeCoeff2TE)
    @test isapprox(reflModeCoeffTE, reflModeCoeffTEbenchmark, rtol=1e-2)
    @test isapprox(transModeCoeffTE, transModeCoeffTEbenchmark, rtol=1e-2)


    # Compute reflected and transmitted fields ("eref" in benchmark)
    EᵣxyTE = Wᵣ*reflModeCoeffTE
    EₜxyTE = Wₜ*transModeCoeffTE
    EᵣxyTEbenchmark = 
        [-0.0094 - 0.0369im; 0.0000 + 0.0000im;    0.0092 + 0.0369im;    0.0012 - 0.0042im;    0.0000 + 0.0000im;    -0.0012 + 0.0042im;    0.0374 + 0.0158im;    0.0000 - 0.0000im;    -0.0373 - 0.0160im;    0.0432 + 0.0301im;    -0.0265 - 0.0191im;    0.0430 + 0.0304im;    -0.1028 + 0.0241im;    -0.2706 - 0.0715im;    -0.1029 + 0.0235im;    0.0378 - 0.0313im;    -0.0199 + 0.0283im; 0.0380 - 0.0310im]
    EₜxyTEbenchmark = 
        [0.0268 + 0.0588im;0.0000 + 0.0000im;-0.0264 - 0.0589im;-0.0134 + 0.0036im;0.0000 + 0.0000im;0.0134 - 0.0035im;-0.0569 - 0.0036im;-0.0000 + 0.0000im;0.0569 + 0.0039im;-0.0263 - 0.0824im;0.0699 + 0.1017im;-0.0258 - 0.0825im;0.1382 - 0.0589im;0.5145 + 0.2249im;0.1386 - 0.0581im;-0.0896 - 0.0165im;0.0809 - 0.0620im;-0.0895 - 0.0171im]
    @test isapprox(EᵣxyTE,EᵣxyTEbenchmark, rtol=1e-2)
    @test isapprox(EₜxyTE,EₜxyTEbenchmark, rtol=1e-2)


    # simpler
    EᵣxyTE = modeCoeff2Field(reflModeCoeffTE, Wᵣ)
    EₜxyTE = modeCoeff2Field(transModeCoeffTE, Wₜ)
    @test isapprox(EᵣxyTE,EᵣxyTEbenchmark, rtol=1e-2)
    @test isapprox(EₜxyTE,EₜxyTEbenchmark, rtol=1e-2)



    # Compute longitudinal field components
    # Unlike reference, these are NOT normalized by default
    KzRefl, KzTrans = calcKzReflTrans(kVectorSet, layerReflection, layerTransmission, matCol)
    KzReflNorm = KzRefl/kVectorSet.k₀  # Normalized z-component
    KzTransNorm = KzTrans/kVectorSet.k₀



    EᵣxTE = EᵣxyTE[1:numHarmonics(kVectorSet)]
    EᵣyTE = EᵣxyTE[(numHarmonics(kVectorSet)+1):(2*numHarmonics(kVectorSet))]
    EᵣzTE = inv(kzᵣ)*(kVectorSet.Kx*EᵣxTE + kVectorSet.Ky*EᵣyTE)

    EₜxTE = EₜxyTE[1:numHarmonics(kVectorSet)]
    EₜyTE = EₜxyTE[(numHarmonics(kVectorSet)+1):(2*numHarmonics(kVectorSet))]
    EₜzTE = inv(kzₜ)*(kVectorSet.Kx*EₜxTE + kVectorSet.Ky*EₜyTE)
    
    @test isapprox(EᵣxTE, EᵣxTEbenchmark, rtol=1e-2)
    @test isapprox(EᵣyTE, EᵣyTEbenchmark, rtol=1e-2)
    @test isapprox(EᵣzTE, EᵣzTEbenchmark, rtol=1e-2)
    @test isapprox(EₜxTE, EₜxTEbenchmark, rtol=1e-2)
    @test isapprox(EₜyTE, EₜyTEbenchmark, rtol=1e-2)
    @test isapprox(EₜzTE, EₜzTEbenchmark, rtol=1e-2)
    


    # EᵣTE = calcFieldsxyz(EᵣxyTE, KzRefl, kVectorSet)
    # EₜTE = calcFieldsxyz(EₜxyTE, KzRefl, kVectorSet)
    EᵣxTE, EᵣyTE, EᵣzTE = calcFieldsxyz(EᵣxyTE, kzᵣ, kVectorSet)
    EₜxTE, EₜyTE, EₜzTE = calcFieldsxyz(EₜxyTE, kzₜ, kVectorSet)
    @test isapprox(EᵣxTE, EᵣxTEbenchmark, rtol=1e-2)
    @test isapprox(EᵣyTE, EᵣyTEbenchmark, rtol=1e-2)
    @test isapprox(EᵣzTE, EᵣzTEbenchmark, rtol=1e-2)
    @test isapprox(EₜxTE, EₜxTEbenchmark, rtol=1e-2)
    @test isapprox(EₜyTE, EₜyTEbenchmark, rtol=1e-2)
    @test isapprox(EₜzTE, EₜzTEbenchmark, rtol=1e-2)



    # Step 12: Diffraction efficiencies
    # Reflected power
    Eᵣ²TE = calcE²(EᵣxTE, EᵣyTE, EᵣzTE)
    reflMat = getMaterial(matCol, layerReflection.backgroundMaterialName)
    μᵣ = calc_μ(reflMat,λ₀)
    nᵣ = calc_n(reflMat,λ₀)
    
    # Only works for single order of incidence.  equation below assumes unit amplitude source.  The kz is used to account for AOI in power
    # RTE = real(k[Z] / μᵣ)/real(k[Z]/μᵣ)*Eᵣ²TE
    kzᵣᵢ = calckzᵣ(kVectorSet, layerReflection, matCol, λ₀)
    Eᵣ²kzᵣ = -kzᵣᵢ.*Eᵣ²TE
    
    EsourcexTE = sourceFields1TE[1:numHarmonics(kVectorSet)]
    EsourceyTE = sourceFields1TE[ (numHarmonics(kVectorSet)+1):(2*numHarmonics(kVectorSet))]
    EsourcezTE = inv(kzᵣ)*(kVectorSet.Kx*EsourcexTE + kVectorSet.Ky*EsourceyTE)
    
    # EsourcexTE, EsourceyTE, EsourceyTE = fieldsXY2fieldsXYZ(sourceFields1TE, kVectorSet)
    
    Esource²TE = calcE²(EsourcexTE, EsourceyTE, EsourcezTE)
    # EsourcexTE, EsourceyTE, EsourcezTE = calcFieldsxyz(EᵣxyTE, kzᵣ, kVectorSet)
    Esource²kz = -kzᵣᵢ.*Esource²TE
    powerInput = sum(real(Esource²kz))
    RTE = real(Eᵣ²kzᵣ) / powerInput
    RTEbenchmark = 
        [0, 0.0032, 0,
        0.0066, 0.0783, 0.0066,
        0, 0.0032, 0]
    @test isapprox(RTE, RTEbenchmark, rtol=1e-2)
    RtotalTE = sum(RTE)
    @test isapprox(RtotalTE, 0.098299, rtol=1e-2)
    

    Eₜ²TE = calcE²(EₜxTE, EₜyTE, EₜzTE)

    transMat = getMaterial(matCol, layerTransmission.backgroundMaterialName)
    μₜ = calc_μ(transMat,λ₀)
    nₜ = calc_n(transMat,λ₀)
    # Only works for single order of incidence.  equation below assumes unit amplitude source.
    # TTE = real(k[Z]*nₜ / μₜ)/real(k[Z]*nᵣ/μₜ)* Eₜ²TE
    TTEbenchmark = 
        [0.0206, 0.0360, 0.0208,
        0.0447, 0.6689, 0.0447,
        0.0206, 0.0360, 0.0208]
        
    kzₜᵢ = calckzₜ(kVectorSet, layerTransmission, matCol, λ₀)
    Eₜ²kzₜ = kzₜᵢ.*Eₜ²TE
    TTE = real(Eₜ²kzₜ) / powerInput
    @test isapprox(TTE, TTEbenchmark, rtol=1e-1)
    TtotalTE = sum(TTE)
    @test isapprox(TtotalTE, 0.8985, rtol=1e-2)

    

    # Straight from source fields to scattered fields
    # Eᵣ, Eₜ = calcScatteredFields(sourceFields1TE, sourceFields2TE, Sglobal, Wᵣ, Wₜ)
    # @test isapprox(Eᵣ,[0;-0.3156-0.0437im], rtol=1e-2)
    # @test isapprox(Eₜ,[0;0.5912+0.2721im], rtol=1e-2)

    return



# TODO MAKE A SIMULATION STRUCT THAT HAS TYPICAL REQUIRED THINGS, KVECTORSET, GVECTORSET, MATDICT

end


end
