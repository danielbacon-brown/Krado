@testset "Material" begin
    material1 = Material(ConstantPermittivity(4), ConstantPermeability(2))
    material2 = Material(ConstantPermittivity(3))
    λ₀ = 1.0 * μm
    wavenumber = WavenumberByλ₀(λ₀)
    @test calc_ϵ(material1, wavenumber) == 4
    @test calc_μ(material1, wavenumber) == 2
    @test calc_ϵ(material2, wavenumber) == 3
    @test calc_μ(material2, wavenumber) == 1
end;
@testset "Function Material" begin
    permittivityFunc(wavenumber) = 2 + (1*μm)^2/getλ₀(wavenumber)^2
    @test permittivityFunc(WavenumberByλ₀(1.0*μm)) ≈ 3.0
    permittivityModel = FunctionPermittivity(permittivityFunc)
    @test calc_ϵ(permittivityModel, WavenumberByλ₀(1.0*μm)) ≈ 3
    material1 = Material( permittivityModel )
    @test calc_ϵ(material1, WavenumberByλ₀(1.0*μm)) ≈ 3
    @test calc_ϵ(material1, WavenumberByλ₀(2.0*μm)) ≈ 2.25
end;
@testset "Spatial Function Material" begin
    permittivityFunc(wavenumber, position) = 1 + sin(position[X]) + cos(position[Y])
    @test permittivityFunc(WavenumberByλ₀(1.0*μm), _2VectorFloat(1/2*pi,2*pi)) ≈ 3.0
    permittivityModel = SpatialFunctionPermittivity(permittivityFunc)
    @test calc_ϵ(permittivityModel, WavenumberByλ₀(1.0*μm), _2VectorFloat(0,0)) ≈ 2.0
    material2 = Material( permittivityModel )
    @test calc_ϵ(material2, WavenumberByλ₀(1.0*μm), _2VectorFloat(1/2*pi, 2*pi)) ≈ 3.0
    # @test calc_ϵ(material1, 2.0 * μm) ≈ 2.25
end;

@testset "Wavenumber limits" begin
    wavenumberRange1 = (WavenumberByλ₀(0.5*μm), WavenumberByλ₀(2*μm) )
    wavenumberRange2 = (WavenumberByλ₀(3*μm), WavenumberByλ₀(1*μm) )
    material1 = Material(ConstantPermittivity(4), ConstantPermeability(2); wavenumberRange = wavenumberRange1)
    material2 = Material(ConstantPermittivity(3); wavenumberRange=wavenumberRange2)
    # λ₀ = 1.0 * μm
    wavenumber = WavenumberByλ₀(1*μm)
    @test calc_ϵ(material1, WavenumberByλ₀(1*μm)) == 4
    @test calc_μ(material1, WavenumberByλ₀(1*μm)) == 2
    @test_throws DomainError calc_ϵ(material1, WavenumberByλ₀(3*μm))
    @test_throws DomainError calc_μ(material1, WavenumberByλ₀(3*μm))
    @test calc_ϵ(material2, WavenumberByλ₀(3*μm)) == 3
    @test calc_μ(material2, WavenumberByλ₀(3*μm)) == 1
    @test_throws DomainError calc_ϵ(material2, WavenumberByλ₀(0.5*μm))
    @test_throws DomainError calc_μ(material2, WavenumberByλ₀(0.5*μm))
end;
