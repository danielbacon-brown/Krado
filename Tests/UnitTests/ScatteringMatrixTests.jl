@testset "Scattering matrices" begin

    Sglobal = initializeGlobalScatteringMatrix(Float64,1)
    @test Sglobal.matrix ≈    [0 0 1 0;
                        0 0 0 1;
                        1 0 0 0;
                        0 1 0 0]

end;
