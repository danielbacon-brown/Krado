@testset "KVectorSet" begin
    # Test calculation of matrices:
    k₀ = 5
    wavenumber = WavenumberByk₀(k₀)
    k₁ = [ _2VectorFloat(1,0), _2VectorFloat(1,2) ]
    kVectorSet = KVectorSet(wavenumber, k₁)
    @test kVectorSet.Kx ≈ [1 0;
                            0 1]

    # Integration:
    # Defining g-vectors:
    lattice = Lattice(2*π)
    ϖᵢ = [ [0,0], [1,0], [2,0] ]
    harmonicsSet = HarmonicsSet(ϖᵢ)
    Gvectors = GvectorSet(harmonicsSet, lattice)

    # Defining incident light:
    k₀ = 3.0
    wavenumber = WavenumberByk₀(k₀)
    k = _2VectorFloat(1.7,0)
    ϖ = _2VectorInt(1,0)
    kVectorSet = createKVectorSet(wavenumber, k, ϖ, Gvectors)
    @test kVectorSet.kᵢ ≈ [[0.7, 0.0], [1.7, 0.0], [2.7, 0.0]]

end
