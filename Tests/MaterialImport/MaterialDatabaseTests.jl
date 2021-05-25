module MaterialDatabaseTestModule

using Test

include("../../src/IncludeKrado.jl")

@testset "MaterialDatabase" begin

matCol = MaterialCollection()

materialDatabasePath = "Tests/MaterialImport/MaterialDatabasePiece/"

include( "UserMaterialImportFunctions.jl")

# Import a tabulated nk text file
path = "Tests/MaterialImport/Ag_JohnsonAndChristy1972.txt"
name = "Ag_JandC"
addMaterial!(matCol, name, importλnkTextMaterial(path; scale=μm, skipRows=1, delimiter="\t"))
wavenumber = WavenumberByλ₀(0.2119*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.2, rtol=1e-3)
@test isapprox(imag(n), 1.325, rtol=1e-3)


# Import YAML files
# tabulated nk
path = materialDatabasePath * "data/main/Ag/Babar.yml"
name = "Ag"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(1.033*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 0.07286, rtol=1e-3)
@test isapprox(imag(n), 7.412, rtol=1e-3)
ϵ = calc_ϵ( getMaterial(matCol,name), wavenumber )
@test isapprox(real(ϵ), -54.932, rtol=1e-3)
@test isapprox(imag(ϵ), 1.0801, rtol=1e-3)
ϵ, μ = calc_ϵμ( getMaterial(matCol,name), wavenumber )
@test isapprox(real(ϵ), -54.932, rtol=1e-3)
@test isapprox(imag(ϵ), 1.0801, rtol=1e-3)
@test isapprox(real(μ), 1, rtol=1e-3)
@test isapprox(imag(μ), 0, rtol=1e-3)
@test isapprox( convert_n2ϵ(n), ϵ, rtol=1e-3)
@test isapprox( convert_ϵ2n(ϵ), n, rtol=1e-3)

@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )

# tabulated n (no k)
path = materialDatabasePath * "data/main/Al2O3/Boidin.yml"
name = "Al2O3_Boidin"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(0.5876*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.6798, rtol=1e-3)
@test isapprox(imag(n), 0, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )

# tabulated n
path = materialDatabasePath * "data/glass/lzos/CTK8.yml"
name = "CTK8"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(0.5876*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.7031, rtol=1e-3)
@test isapprox(imag(n), 0, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )


# formula 1
path = materialDatabasePath * "data/main/Al2O3/Malitson.yml"
name = "Al2O3_Malitson"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(1*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.7557, rtol=1e-3)
@test isapprox(imag(n), 0, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )


# formula 2
path = materialDatabasePath * "data/glass/Schott/N-BK7.yml"
name = "BK7"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(0.5876*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.5168, rtol=1e-3)
@test isapprox(imag(n), 9.7525e-9, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )

# formula 3
path = materialDatabasePath * "data/glass/cdgm/BAF2.yml"
name = "BAF2"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(0.5876*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.5697, rtol=1e-3)
@test isapprox(imag(n), 1.4049e-8, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )

# formula 4
path = materialDatabasePath * "data/organic/CH4N2O - urea/Rosker-o.yml"
name = "Urea-no"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(0.5876*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.4906, rtol=1e-3)
@test isapprox(imag(n), 0, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )

# formula 5
path = materialDatabasePath * "data/organic/C2H6 - ethane/Loria.yml"
name = "Ethane"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(0.523*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.00075794, rtol=1e-5)
@test isapprox(imag(n), 0, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )

# formula 6
path = materialDatabasePath * "data/main/N2/Griesmann.yml"
name = "N2_Griesmann"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(0.15*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.00039622, rtol=1e-6)
@test isapprox(imag(n), 0, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )

# formula 7
path = materialDatabasePath * "data/main/Si/Edwards.yml"
name = "Si_Edwards"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(2.4373*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 3.4434, rtol=1e-3)
@test isapprox(imag(n), 0, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )


# formula 8
path = materialDatabasePath * "data/main/AgBr/Schroter.yml"
name = "AgBr"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(0.5876*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 2.2579, rtol=1e-3)
@test isapprox(imag(n), 0, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )

# formula 9
path = materialDatabasePath * "data/organic/CH4N2O - urea/Rosker-e.yml"
name = "Urea -ne"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(0.5876*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.6065, rtol=1e-3)
@test isapprox(imag(n), 0, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )


# Air
path = materialDatabasePath * "data/other/mixed gases/air/Ciddor.yml"
name = "Air_Ciddor"
addMaterial!(matCol, name, importYAMLmaterial(path))
wavenumber = WavenumberByλ₀(0.5876*μm)
n = calc_n( getMaterial(matCol,name), wavenumber )
@test isapprox(real(n), 1.00027717, rtol=1e-5)
@test isapprox(imag(n), 0, rtol=1e-3)
@test_throws DomainError calc_ϵμ( getMaterial(matCol,name), WavenumberByλ₀(1000*μm)  )


# Load all materials:
# data\main\ZnS\Amotchkina.yml must be changed for this to work
# for (root, dirs, files) in walkdir("MaterialDatabase\\data")
#     for file in files
#         path = joinpath(root, file)
#         println(joinpath(root, file)) # path to files
#         extension = file[findlast(isequal('.'),file):end]
#         if extension == ".yml"
#             addMaterial!(matCol, file, importYAMLmaterial(path))
#         end
#     end
# end



# Test the one-word import of materials
userMaterialImportFunctions = getUserMaterialImportFunctions()
matCol= MaterialCollection()
importUserMaterial!(matCol, userMaterialImportFunctions, "Ag_J&C")
wavenumber = WavenumberByλ₀(0.2119*μm)
n = calc_n( getMaterial(matCol,"Ag_J&C"), wavenumber )
@test isapprox(real(n), 1.2, rtol=1e-3)
@test isapprox(imag(n), 1.325, rtol=1e-3)


# Test getMaterialsUsed in stack
substrateLayer = SemiInfiniteLayerDefinition("Al2O3")
superstrateLayer = SemiInfiniteLayerDefinition("Air")
layer1 = UniformLayerDefinition(50*nm, "Ag_J&C")
layer2solids = [Solid(Rectangle([0,0],[0.518*μm,0.324*μm]), "Al2O3")]
layer2 = PatternedLayerDefinition([100,100], 100*nm, LayerPattern(layer2solids, "Air"))
layer3solids = [Solid(Circle([0.1*μm,0],0.15*μm), "Ag_J&C")]
layer3 = PatternedLayerDefinition([100,100], 100*nm, LayerPattern(layer3solids, "Al2O3"))
layer4solids = [Solid(Polygon( [[0,0],[0.1,0],[-0.05,0.05],[0,-0.1]]*μm; offset=[0.1*μm,0]), "Ag_J&C")]
layer4 = PatternedLayerDefinition([100,100], 100*nm, LayerPattern(layer4solids, "Air"))
layerStack = [substrateLayer, layer1, layer2, layer3, layer4, superstrateLayer]

materialsUsed = getMaterialsUsed(layerStack)
@test materialsUsed == Set{String}(["Air", "Al2O3", "Ag_J&C"])

# Test importing of getMaterialsUsed list
matCol = MaterialCollection()
userMaterialImportFunctions = getUserMaterialImportFunctions()
importUserMaterial!(matCol,  userMaterialImportFunctions, getMaterialsUsed(layerStack))
n = calc_n( getMaterial(matCol,"Ag_J&C"), wavenumber )
@test isapprox(real(n), 1.2, rtol=1e-3)
@test isapprox(imag(n), 1.325, rtol=1e-3)

# Test one-line importing of getMaterialsUsed list
matCol = importUserMaterial( getUserMaterialImportFunctions(), getMaterialsUsed(layerStack))
@test contains(matCol,"Vacuum")  == false  # Not included because it was not in the simulation
n = calc_n( getMaterial(matCol,"Ag_J&C"), wavenumber )
@test isapprox(real(n), 1.2, rtol=1e-3)
@test isapprox(imag(n), 1.325, rtol=1e-3)


end;


end
