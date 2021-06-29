
function getUserMaterialImportFunctions()
    userMaterialImportFunctions = Dict{Any,Any}()

    userMaterialImportFunctions["Vacuum"] = ()->Material(ConstantPermittivity(1))
    userMaterialImportFunctions["Air"] = ()->Material(ConstantPermittivity(1))
    userMaterialImportFunctions["Si Endura"] = ()->importλnkTextMaterial("C:/UserMaterials/Si_Endura.txt"; scale=μm, skipRows=1, delimiter=" ")
    # userMaterialImportFunctions["Ag_J&C"] = ()->importλnkTextMaterial("Tests/MaterialImport/Ag_JohnsonAndChristy1972.txt"; scale=μm, skipRows=1, delimiter="\t")

    return userMaterialImportFunctions
end


function getUserPlottingParameters()
    # plottingParameters = Dict{Any,PlottingParameters}()
    plottingParameters = PlottingParameterCollection()
    addPlottingParameter!( plottingParameters, "Vacuum", PlottingParameters(;
            color = [0.99, 0.99, 0.99],
            alpha = 0,
            shade = false,
            lineWidth = 0,
            lineStyle = "None",
        ) )
    addPlottingParameter!( plottingParameters, "Air", PlottingParameters(;
            color = [0.95, 0.95, 0.95],
            alpha = 0,
            shade = false,
            lineWidth = 0,
            lineStyle = "None",
        ) )
    addPlottingParameter!( plottingParameters, "Si Endura", PlottingParameters(;
            color = [0.7, 0.7, 0.4],
            alpha = 0.1,
            shade = true,
            lineStyle = "-",
            lineColor = "r",
        ) )
    addPlottingParameter!( plottingParameters, "Ag", PlottingParameters(;
            color = [0.4, 0.4, 0.4],
            alpha = 0.25,
            shade = true,
            lineStyle = "-",
            lineColor = "b",
        ) )
    addPlottingParameter!( plottingParameters, "Al2O3", PlottingParameters(;
            color = [0.2, 0.2, 0.7],
            alpha = 0.25,
            shade = true,
            lineStyle = "-",
            lineColor = [0.2, 0.7, 0.2],
        ) )
    return plottingParameters
end
