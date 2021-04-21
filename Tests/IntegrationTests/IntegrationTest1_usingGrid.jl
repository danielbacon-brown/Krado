module IntegrationTest1
# Using RCWA-Benchmark-Data-1x1 by Raymond Rumpf from empossible.net

using Test
# using LinearAlgebra
# using StaticArrays
# import Statistics

println()
println()

@testset "IntegrationTest1_usingGrid" begin

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
    # Using 1-order only:
    M = 0
    N = 0
    harmonicsDef = HarmonicsTruncationByRectangle(M,N)
    harmonicsSet = calcHarmonicsSet(harmonicsDef)
    @test harmonicsSet.mnᵢ == [_2VectorInt(0, 0)]
    @test harmonicsSet.indᵢ_mn == Dict(_2VectorInt(0, 0)=>1)

    #K-vectors and G-vectors
    Gvectors = GvectorSet(harmonicsSet,lattice)


    # Define layers.  Layer 1 contains triangle
    layer1_d = 0.5 * cm  # thickness of layer 1
    layer2_d = 0.3 * cm

    # Define unpatterned layer:
    layer2 = UniformLayerDefinition(layer2_d, "device_material")


    # Define shape: equilateral triangle pointed toward +y.

    # simple style:
    triangle_length = 0.8 * Ly
    # point1 = [triangle_length/2, -(triangle_length/2) / tand(60)]
    # point2 = [-triangle_length/2, -(triangle_length/2) / tand(60)]
    # point3 = [0, (triangle_length/2) / sind(60)]

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

    # Define using a triangle, the preferred method, but results only approximate the benchmark.
    # point1 = [triangle_length/2, -(triangle_length/2) / tand(60)]
    # point2 = [-triangle_length/2, -(triangle_length/2) / tand(60)]
    # point3 = [0, (triangle_length/2) / sind(60)]
    # triangle_shape = Polygon([Lx/2,Ly/2],[point1, point2, point3])
    # triangle_solid = Solid( triangle_shape, "refl_material" )
    # Define solids calculator
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
    @test isapprox( Cϵᵢⱼ1[1,1], 5.0449, rtol=1e-4)
    @test Cμᵢⱼ1[1,1] ≈ Complex(1,0)

    Cϵᵢⱼ2, Cμᵢⱼ2 = calcConvolutionMatrices( layer2, lattice, Gvectors, matCol, λ₀ )
    @test Cϵᵢⱼ2[1,1] ≈ Complex(6,0)
    @test Cμᵢⱼ2[1,1] ≈ Complex(1,0)

    #Calculate inverse convolution matrices for each layer
    Cϵᵢⱼ⁻¹1 = inv(Cϵᵢⱼ1)
    Cμᵢⱼ⁻¹1 = inv(Cμᵢⱼ1)
    Cϵᵢⱼ⁻¹2 = inv(Cϵᵢⱼ2)
    Cμᵢⱼ⁻¹2 = inv(Cμᵢⱼ2)





    # Wave vector expansion:
    @test isapprox( convert_ϵ2n(calc_ϵ( getMaterial(matCol,"refl_material"), λ₀)), 1.4142, rtol = 1e-3)
    @test isapprox( convert_ϵ2n(calc_ϵ( getMaterial(matCol,"tran_material"), λ₀)), 3, rtol = 1e-3)
    @test isapprox( λ₀2k₀(λ₀), 314.1593, rtol = 1e-3)

    # Create the incident wave-vector
    k₀ = λ₀2k₀(λ₀)
    k = normalIncidenceKvector(k₀)
    @test norm(k) ≈ k₀
    kVectorSet = createKVectorSet(k, Gvectors)

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

    # Verify the harmonic orders (denoted p,q in the lecture materials)
    @test harmonicsSet.mnᵢ == [_2VectorInt(0,0)]


    #STEP 5: Calc eigenmodes for free space
    # same for all layers
    KzNorm = calcKzForUnpatternedEigenmode(kVectorSet)
    @test KzNorm ≈ [1]

    Q = calcQForUnpatternedEigenmode( kVectorSet )
    @test Q == [0 1; -1 0]

    W₀ = calcW₀( numHarmonics(kVectorSet) )
    @test W₀ ≈ [1 0; 0 1]

    λeigenvalues = calcλForUnpatternedEigenMode(KzNorm)
    @test λeigenvalues ≈ [ im*1 0;
                        0 im*1]

    V₀ = calcEigenmodesFromQλ(Q,λeigenvalues)
    @test V₀ ≈ [ 0  im*-1;
                im*1  0]

    V₀ = calcV₀(kVectorSet)
    @test V₀ ≈ [ 0  im*-1;
                im*1  0]

    #STEP 6: Initialize global S
    Sglobal = initializeGlobalS(numHarmonics(harmonicsSet))
    @test Sglobal ≈    [0 0 1 0;
                        0 0 0 1;
                        1 0 0 0;
                        0 1 0 0]


    # STEP 7: MAIN LOOP

    # Layer 1: Patterned
    P₁ = calcPmatrixPatterned(kVectorSet, Cϵᵢⱼ1, Cϵᵢⱼ⁻¹1, Cμᵢⱼ1, Cμᵢⱼ⁻¹1)
    @test P₁ ≈ [0 1;
        -1 0]
    Q₁ = calcQmatrixPatterned(kVectorSet, Cϵᵢⱼ1, Cϵᵢⱼ⁻¹1, Cμᵢⱼ1, Cμᵢⱼ⁻¹1)
    @test isapprox(Q₁, [0 5.0449;
        -5.0449 0], rtol=1e-3)

    Ω²₁ = calcΩ²(P₁,Q₁)
    @test isapprox(Ω²₁, [-5.0449 0;
                        0 -5.0449], rtol=1e-3)

    W₁, λ₁ = calcWᵢλᵢ(Ω²₁)
    @test W₁ ≈ [1 0;
                0 1]
    @test isapprox( λ₁, [2.2461im 0;
                0 2.2461im], rtol = 1e-3)

    V₁ = calcEigenmodesFromQλ(Q₁,λ₁)
    @test isapprox( V₁, [0 -2.246im;
                    2.246im 0], rtol = 1e-3)

    A₁, B₁ = calcAB(W₁,W₀,V₁,V₀)
    @test isapprox(A₁, [1.4452 0;
                        0 1.4452], rtol=1e-3)
    @test isapprox(B₁, [0.5548 0;
                        0 0.5548], rtol=1e-3)
    X₁ = calcX(λ₁, kVectorSet.k₀, layer1.thickness)
    @test isapprox(X₁, [-0.9262+0.3770im 0;
                        0 -0.9262+0.3770im], rtol=1e-3)

    _1, _2 = getQuadrantSlices(numHarmonics(kVectorSet))

    S₁ = calcScatteringMatrix_ABX(A₁,B₁,X₁)
    @test isapprox(S₁[_1,_1], [(-0.1544 - 0.2819im) (0);
                                (0) (-0.1544 - 0.2819im)], rtol=1e-3)
    @test isapprox(S₁[_1,_2], [(-0.8305 + 0.4549im) (0);
                                (0) (-0.8305 + 0.4549im)], rtol=1e-3)
    @test isapprox(S₁[_2,_1], [(-0.8305 + 0.4549im) (0);
                                (0) (-0.8305 + 0.4549im)], rtol=1e-3)
    @test isapprox(S₁[_2,_2], [(-0.1544 - 0.2819im) (0);
                                (0) (-0.1544 - 0.2819im)], rtol=1e-3)


    # sugary:
    S₁ = calcScatteringMatrix(layer1, matCol, kVectorSet, Gvectors, lattice )
    @test isapprox(S₁[_1,_1], [(-0.1544 - 0.2819im) (0);
                                (0) (-0.1544 - 0.2819im)], rtol=1e-3)
    @test isapprox(S₁[_1,_2], [(-0.8305 + 0.4549im) (0);
                                (0) (-0.8305 + 0.4549im)], rtol=1e-3)
    @test isapprox(S₁[_2,_1], [(-0.8305 + 0.4549im) (0);
                                (0) (-0.8305 + 0.4549im)], rtol=1e-3)
    @test isapprox(S₁[_2,_2], [(-0.1544 - 0.2819im) (0);
                                (0) (-0.1544 - 0.2819im)], rtol=1e-3)


    # Layer 2: Unpatterned
    mat₂str = layer2.backgroundMaterialName
    # mat₂ = matDict[mat₂str]
    mat₂ = getMaterial(matCol,mat₂str)
    ϵ₂, μ₂ = calc_ϵμ(mat₂,λ₀)
    P₂ = calcPmatrixUnpatterned(kVectorSet, ϵ₂, μ₂ )
    @test P₂ ≈ [ 0 1;
                -1 0]
    Q₂ = calcQmatrixUnpatterned(P₂, ϵ₂, μ₂)
    @test Q₂ ≈ [ 0 6;
                -6 0]

    # Sugary version.  Method specific to uniform layer
    # P₂, Q₂ = calcPQmatrix(layer2, kVectorSet, matDict, λ₀)
    # P₂, Q₂ = calcPQmatrix(layer2, kVectorSet, matDict)
    P₂, Q₂ = calcPQmatrix(layer2, kVectorSet, matCol)
    @test P₂ ≈ [ 0 1;
                -1 0]
    @test Q₂ ≈ [ 0 6;
                -6 0]

    # calcΩ² same for each layer
    Ω²₂ = calcΩ²(P₂,Q₂)
    @test Ω²₂ ≈ [-6  0;
                  0 -6]

    W₂, λ₂ = calcWᵢλᵢ(Ω²₂)
    @test W₂ ≈ [1 0;
                0 1]
    @test isapprox( λ₂, [2.4495im 0;
                0 2.4495im], rtol = 1e-3)

    V₂ = calcEigenmodesFromQλ(Q₂,λ₂)
    @test isapprox( V₂, [0 -2.4495im;
                2.4495im 0], rtol = 1e-3)

    # sugary:
    V₂ = calcEigenmodesForUniformLayer(kVectorSet, layer2, matCol)
    @test isapprox( V₂, [0 -2.4495im;
                2.4495im 0], rtol = 1e-3)


    # Common components of scattering matrix
    A₂ = calcA(W₂, W₀, V₂, V₀)
    @test isapprox(A₂, [1.4082 0;
                0 1.4082], rtol=1e-3)
    B₂ = calcB(W₂, W₀, V₂, V₀)
    @test isapprox(B₂, [0.5918 0;
                0 0.5918], rtol=1e-3)

    X₂ = calcX(λ₂, kVectorSet.k₀, layer2.thickness)
    @test isapprox( X₂, [(-0.6727-0.7400im) (0);
                        (0) (-0.6727-0.7400im)], rtol=1e-3 )

    S₂ = calcScatteringMatrix_ABX(A₂, B₂, X₂)

    @test isapprox(S₂[_1,_1], [(-0.5085 + 0.3235im) (0);
                                (0)   (-0.5085 + 0.3235im)], rtol=1e-3)
    @test isapprox(S₂[_1,_2], [(-0.4283 - 0.6733im) (0);
                                (0)   (-0.4283 - 0.6783im)], rtol=1e-2)
    @test isapprox(S₂[_2,_1], [(-0.4283 - 0.6733im) (0);
                                (0)   (-0.4283 - 0.6733im)], rtol=1e-2)
    @test isapprox(S₂[_2,_2], [(-0.5085 + 0.3235im) (0);
                                (0)   (-0.5085 + 0.3235im)], rtol=1e-3)


    # S₂ = calcScatteringMatrixUniform(layer2, matDict, kVectorSet)
    S₂ = calcScatteringMatrix(layer2, matCol, kVectorSet)
    @test isapprox(S₂[_1,_1], [(-0.5085 + 0.3235im) (0);
                                (0)   (-0.5085 + 0.3235im)], rtol=1e-3)
    @test isapprox(S₂[_1,_2], [(-0.4283 - 0.6733im) (0);
                                (0)   (-0.4283 - 0.6783im)], rtol=1e-2)
    @test isapprox(S₂[_2,_1], [(-0.4283 - 0.6733im) (0);
                                (0)   (-0.4283 - 0.6733im)], rtol=1e-2)
    @test isapprox(S₂[_2,_2], [(-0.5085 + 0.3235im) (0);
                                (0)   (-0.5085 + 0.3235im)], rtol=1e-3)


    # STEP 8: Reflection side scattering matrix
    Pᵣ, Qᵣ = calcPQmatrix(layerReflection, kVectorSet, matCol)
    @test Qᵣ ≈ [0 2;
                -2 0]
    Ω²ᵣ = calcΩ²(Pᵣ,Qᵣ)

    Wᵣ, λᵣ = calcWᵢλᵢ(Ω²ᵣ)
    @test isapprox(λᵣ, [1.4142im 0;
                0 1.4142im], rtol=1e-3)
    @test Wᵣ ≈ [1 0;
                0 1]

    Vᵣ = calcEigenmodesFromQλ(Qᵣ,λᵣ)
    @test isapprox(Vᵣ, [0 -1.4142im;
                1.4142im 0], rtol=1e-3)

    Aᵣ = calcA_SemiInfinite(Wᵣ, W₀, Vᵣ, V₀)
    @test isapprox(Aᵣ, [2.4142 0;
                0 2.4142], rtol=1e-3)
    Bᵣ = calcB_SemiInfinite(Wᵣ, W₀, Vᵣ, V₀)
    @test isapprox(Bᵣ, [-0.4142 0;
                0 -0.4142], rtol=1e-3)
    # #sugary
    Aᵣ, Bᵣ = calcAB_SemiInfinite(Wᵣ, W₀, Vᵣ, V₀)
    @test isapprox(Aᵣ, [2.4142 0;
                0 2.4142], rtol=1e-3)
    @test isapprox(Bᵣ, [-0.4142 0;
                0 -0.4142], rtol=1e-3)

    Sᵣ = calcScatteringMatrixReflection_AB(Aᵣ,Bᵣ)
    @test isapprox(Sᵣ[_1,_1], [(0.1716) (0);
                                (0)   (0.1716)], rtol=1e-3)
    @test isapprox(Sᵣ[_1,_2], [(0.8284) (0);
                                (0)   (0.8284)], rtol=1e-2)
    @test isapprox(Sᵣ[_2,_1], [(1.1716) (0);
                                (0)   (1.1716)], rtol=1e-2)
    @test isapprox(Sᵣ[_2,_2], [(-0.1716) (0);
                                (0)   (-0.1716)], rtol=1e-3)

    Sᵣ = calcScatteringMatrixReflection(layerReflection, matCol, kVectorSet)
    @test isapprox(Sᵣ[_1,_1], [(0.1716) (0);
                                (0)   (0.1716)], rtol=1e-3)
    @test isapprox(Sᵣ[_1,_2], [(0.8284) (0);
                                (0)   (0.8284)], rtol=1e-2)
    @test isapprox(Sᵣ[_2,_1], [(1.1716) (0);
                                (0)   (1.1716)], rtol=1e-2)
    @test isapprox(Sᵣ[_2,_2], [(-0.1716) (0);
                                (0)   (-0.1716)], rtol=1e-3)


    # STEP 9: Transmission side scattering matrix
    Pₜ, Qₜ = calcPQmatrix(layerTransmission, kVectorSet, matCol)
    @test Qₜ ≈ [0 9;
                -9 0]
    Ω²ₜ = calcΩ²(Pₜ,Qₜ)

    Wₜ, λₜ = calcWᵢλᵢ(Ω²ₜ)
    @test isapprox(λₜ, [3im 0;
                0 3im], rtol=1e-3)
    @test Wₜ ≈ [1 0;
                0 1]

    Vₜ = calcEigenmodesFromQλ(Qₜ,λₜ)
    @test isapprox(Vₜ, [0 -3im;
                3im 0], rtol=1e-3)

    Aₜ = calcA_SemiInfinite(Wₜ, W₀, Vₜ, V₀)
    @test isapprox(Aₜ, [4 0
                0 4], rtol=1e-3)
    Bₜ = calcB_SemiInfinite(Wₜ, W₀, Vₜ, V₀)
    @test isapprox(Bₜ, [-2 0
                0 -2], rtol=1e-3)

    Sₜ = calcScatteringMatrixTransmission_AB(Aₜ,Bₜ)
    @test isapprox(Sₜ[_1,_1], [(-0.5) (0);
                                (0)   (-0.5)], rtol=1e-3)
    @test isapprox(Sₜ[_1,_2], [(1.5) (0);
                                (0)   (1.5)], rtol=1e-2)
    @test isapprox(Sₜ[_2,_1], [(0.5) (0);
                                (0)   (0.5)], rtol=1e-2)
    @test isapprox(Sₜ[_2,_2], [(0.5) (0);
                                (0)   (0.5)], rtol=1e-3)

    Sₜ = calcScatteringMatrixTransmission(layerTransmission, matCol, kVectorSet)
    @test isapprox(Sₜ[_1,_1], [(-0.5) (0);
                                (0)   (-0.5)], rtol=1e-3)
    @test isapprox(Sₜ[_1,_2], [(1.5) (0);
                                (0)   (1.5)], rtol=1e-2)
    @test isapprox(Sₜ[_2,_1], [(0.5) (0);
                                (0)   (0.5)], rtol=1e-2)
    @test isapprox(Sₜ[_2,_2], [(0.5) (0);
                                (0)   (0.5)], rtol=1e-3)

    # STEP 10: Global scattering matrix:
    Sdevice = S₁⊗S₂
    Sglobal = Sᵣ⊗Sdevice⊗Sₜ

    @test isapprox(Sglobal[_1,_1], [(-0.3156-0.0437im) (0);
                                    (0) (-0.3156-0.0437im)], rtol=1e-3)
    @test isapprox(Sglobal[_1,_2], [(1.2541+0.5773im) (0);
                                    (0) (1.2541+0.5773im)], rtol=1e-3)
    @test isapprox(Sglobal[_2,_1], [(0.5912+0.2721im) (0);
                                    (0) (0.5912+0.2721im)], rtol=1e-3)
    @test isapprox(Sglobal[_2,_2], [(0.2384+0.2113im) (0);
                                    (0) (0.2384+0.2113im)], rtol=1e-3)

    # Put it in terms of a device stack:
    layerStack = [layerReflection, layer1, layer2, layerTransmission]
    # Sglobal = calcScatteringMatrix(simulation, layerStack)
    Sglobal = calcScatteringMatrix(layerStack, matCol, kVectorSet, Gvectors, lattice)

    @test isapprox(Sglobal[_1,_1], [(-0.3156-0.0437im) (0);
                                    (0) (-0.3156-0.0437im)], rtol=1e-3)
    @test isapprox(Sglobal[_1,_2], [(1.2541+0.5773im) (0);
                                    (0) (1.2541+0.5773im)], rtol=1e-3)
    @test isapprox(Sglobal[_2,_1], [(0.5912+0.2721im) (0);
                                    (0) (0.5912+0.2721im)], rtol=1e-3)
    @test isapprox(Sglobal[_2,_2], [(0.2384+0.2113im) (0);
                                    (0) (0.2384+0.2113im)], rtol=1e-3)


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
    @test sourceFields1TE ≈ [0;1] # top is x-direction and bottom is y-direction
    sourceFields1TM, sourceFields2TM = calcSourceFields(boundaryConditionsTM, harmonicsSet, kVectorSet)
    @test sourceFields1TM ≈ [1;0] # top is x-direction and bottom is y-direction


    # sourceModeCoeff = inv(Wᵣ)*sourceFieldsTE # Wref^-1 * sourceFieldsTE
    sourceModeCoeff1TE = modeField2Coeff(sourceFields1TE, Wᵣ)
    sourceModeCoeff2TE = modeField2Coeff(sourceFields2TE, Wₜ)

    S₁₁ = Sglobal[_1,_1]
    S₁₂ = Sglobal[_1,_2]
    S₂₁ = Sglobal[_2,_1]
    S₂₂ = Sglobal[_2,_2]

    # Compute transmission and reflection mode coefficients
    reflModeCoeffTE = S₁₁*sourceModeCoeff1TE
    transModeCoeffTE = S₂₁*sourceModeCoeff1TE
    @test isapprox(reflModeCoeffTE, [0;-0.3156-0.0437im], rtol=1e-3)
    @test isapprox(transModeCoeffTE, [0;0.5912+0.2721im], rtol=1e-3)

    # simpler:
    reflModeCoeffTE, transModeCoeffTE = propagateModeCoeff(Sglobal, sourceModeCoeff1TE, sourceModeCoeff2TE)
    @test isapprox(reflModeCoeffTE, [0;-0.3156-0.0437im], rtol=1e-3)
    @test isapprox(transModeCoeffTE, [0;0.5912+0.2721im], rtol=1e-3)


    # Compute reflected and transmitted fields
    EᵣxyTE = Wᵣ*reflModeCoeffTE
    EₜxyTE = Wₜ*transModeCoeffTE
    @test isapprox(EᵣxyTE,[0;-0.3156-0.0437im], rtol=1e-3)
    @test isapprox(EₜxyTE,[0;0.5912+0.2721im], rtol=1e-3)

    # simpler
    EᵣxyTE = modeCoeff2Field(reflModeCoeffTE, Wᵣ)
    EₜxyTE = modeCoeff2Field(transModeCoeffTE, Wₜ)
    @test isapprox(EᵣxyTE,[0;-0.3156-0.0437im], rtol=1e-3)
    @test isapprox(EₜxyTE,[0;0.5912+0.2721im], rtol=1e-3)

    # Compute longitudinal field components
    # Unlike reference, these are NOT normalized by default
    KzRefl, KzTrans = calcKzReflTrans(kVectorSet, layerReflection, layerTransmission, matCol)
    KzReflNorm = KzRefl/kVectorSet.k₀  # Normalized z-component
    KzTransNorm = KzTrans/kVectorSet.k₀
    @test isapprox(KzRefl/kVectorSet.k₀, [-1.4142], rtol=1e-3)
    @test isapprox(KzTrans/kVectorSet.k₀, [3], rtol=1e-3)
    @test isapprox(KzReflNorm, [-1.4142], rtol=1e-3)
    @test isapprox(KzTransNorm, [3], rtol=1e-3)

    EᵣxTE = EᵣxyTE[X,:]
    EᵣyTE = EᵣxyTE[Y,:]
    EᵣzTE = -inv(KzRefl)*(kVectorSet.Kx*EᵣxTE + kVectorSet.Ky*EᵣyTE)
    EₜxTE = EₜxyTE[X,:]
    EₜyTE = EₜxyTE[Y,:]
    EₜzTE = -inv(KzTrans)*(kVectorSet.Kx*EₜxTE + kVectorSet.Ky*EₜyTE)
    @test isapprox(EᵣxTE,[0], rtol=1e-3)
    @test isapprox(EᵣyTE,[-0.3156-0.0437im], rtol=1e-3)
    @test isapprox(EᵣzTE,[0], rtol=1e-3)
    @test isapprox(EₜxTE,[0], rtol=1e-3)
    @test isapprox(EₜyTE,[0.5912+0.2721im], rtol=1e-3)
    @test isapprox(EₜzTE,[0], rtol=1e-3)

    EᵣTE = calcFieldsxyz(EᵣxyTE, KzRefl, kVectorSet)
    EₜTE = calcFieldsxyz(EₜxyTE, KzRefl, kVectorSet)
    @test isapprox(EᵣTE[X],0, rtol=1e-3)
    @test isapprox(EᵣTE[Y],-0.3156-0.0437im, rtol=1e-3)
    @test isapprox(EᵣTE[Z],0, rtol=1e-3)
    @test isapprox(EₜTE[X],0, rtol=1e-3)
    @test isapprox(EₜTE[Y],0.5912+0.2721im, rtol=1e-3)
    @test isapprox(EₜTE[Z],0, rtol=1e-3)

    # Step 12: Diffraction efficiencies
    # Reflected power
    # Eᵣ²TE = abs.(EᵣTE[X,:]).^2 + abs.(EᵣTE[Y,:]).^2 + abs.(EᵣTE[Z,:]).^2
    Eᵣ²TE = calcE²(EᵣTE)

    reflMat = getMaterial(matCol, layerReflection.backgroundMaterialName)
    μᵣ = calc_μ(reflMat,λ₀)
    nᵣ = calc_n(reflMat,λ₀)
    # Only works for single order of incidence.  equation below assumes unit amplitude source.  The kz is used to account for AOI in power
    RTE = real(k[Z] / μᵣ)/real(k[Z]/μᵣ)*Eᵣ²TE
    RtotalTE = sum(RTE)
    @test isapprox(RTE, [0.1015], rtol=1e-3)

    Eₜ²TE = abs.(EₜTE[X,:]).^2 + abs.(EₜTE[Y,:]).^2 + abs.(EₜTE[Z,:]).^2
    # absEₜ²TE = abs.(EₜxTE).^2 + abs.(EₜyTE).^2 + abs.(EₜzTE).^2
    transMat = getMaterial(matCol, layerTransmission.backgroundMaterialName)
    μₜ = calc_μ(transMat,λ₀)
    nₜ = calc_n(transMat,λ₀)
    # Only works for single order of incidence.  equation below assumes unit amplitude source.
    TTE = real(k[Z]*nₜ / μₜ)/real(k[Z]*nᵣ/μₜ)* Eₜ²TE
    TtotalTE = sum(TTE)
    @test isapprox(TTE, [0.8985], rtol=1e-3)


    # Straight from source fields to scattered fields
    Eᵣ, Eₜ = calcScatteredFields(sourceFields1TE, sourceFields2TE, Sglobal, Wᵣ, Wₜ)
    @test isapprox(Eᵣ,[0;-0.3156-0.0437im], rtol=1e-3)
    @test isapprox(Eₜ,[0;0.5912+0.2721im], rtol=1e-3)





# TODO MAKE A SIMULATION STRUCT THAT HAS TYPICAL REQUIRED THINGS, KVECTORSET, GVECTORSET, MATDICT

end


end
