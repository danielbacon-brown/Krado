
function getUserMaterialImportFunctions()
    userMaterialImportFunctions = Dict{Any,Any}()

    userMaterialImportFunctions["Vacuum"] = ()->Material(ConstantPermittivity(1))
    userMaterialImportFunctions["Ag_J&C"] = ()->importλnkTextMaterial("Tests/MaterialImport/Ag_JohnsonAndChristy1972.txt"; scale=μm, skipRows=1, delimiter="\t")

    return userMaterialImportFunctions
end
