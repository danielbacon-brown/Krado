# Material database.
# Refractiveindex.info is a recommended source for optical properties of materials.  This software is designed to work with the database found at: https://refractiveindex.info/about
# All are assumed to have permeability = 1
# If there is one data-section of the YML, and it is for type "formula 1", then parseFormula1 function will be used, and k is set to 0.  If there are two dataType, and the first is for type "formula 1", then parseFormula1 will be used for the value of n, and different function will be used for k.



function importλnkTextMaterial(path::String; scale=μm, skipRows=0, delimiter="\t")

    # Create parallel vectors defining λ₀, n, k
    λ₀Values = Float64[]
    nRealValues = Float64[]
    nImagValues = Float64[]

    open(path) do file
        # Skip header lines
        for iRow in 1:skipRows
            readline(file)
        end

        for line in eachline(file)
            values = split(line, delimiter)
            push!(λ₀Values, parse(Float64,values[1]) )
            push!(nRealValues, parse(Float64,values[2]) )
            push!(nImagValues, parse(Float64,values[3]) )
        end
    end

    # Convert to complex n, and put in SI units.
    nValues = nRealValues .+ nImagValues*1im
    λ₀Values = λ₀Values * scale

    # Sort
    λ₀ValuesSorted, nValuesSorted = sortParallelVectors(λ₀Values, nValues)


    # Create function to interpolate at a given wavenumber
    λnInterpolation = LinearInterpolation(λ₀ValuesSorted, nValuesSorted)
    function interpolateFunction(wavenumber::Wavenumber)::ComplexF64
        return λnInterpolation(getλ₀(wavenumber))
    end
    ϵFunction = x -> convert_n2ϵ(interpolateFunction(x))
    permittivity = FunctionPermittivity(ϵFunction)

    minλ₀ = minimum(λ₀ValuesSorted)
    maxλ₀ = maximum(λ₀ValuesSorted)
    wavenumberRange = (WavenumberByλ₀(minλ₀), WavenumberByλ₀(maxλ₀))

    return Material(permittivity; wavenumberRange = wavenumberRange )
end

function importYAMLmaterial(path::String)


    allData = YAML.load(open(path))
    subData = allData["DATA"]

    # Create a function defining the permittivity at a particular wavenumber.
    if length(subData) == 2  # n and k are separate
        # Get the functions defining n and k from the dictionary based on the "type" label.
        nDataType = subData[1]["type"]
        nFunc, minλ₀n, maxλ₀n = MATERIALDATABASEFUNCTIONS[nDataType](subData[1])
        kDataType = subData[2]["type"]
        kFunc, minλ₀k, maxλ₀k = MATERIALDATABASEFUNCTIONS[kDataType](subData[2])
        ϵFunction = x -> convert_n2ϵ(nFunc(x) + 1im*kFunc(x))
        permittivity = FunctionPermittivity(ϵFunction)
        minλ₀ = max(minλ₀n, minλ₀k)
        maxλ₀ = min(maxλ₀n, maxλ₀k)
    else  # n and k are defined together
        # Get the function defining n and k from the dictionary based on the "type" label.
        dataType = subData[1]["type"]
        nkFunc, minλ₀, maxλ₀ = MATERIALDATABASEFUNCTIONS[dataType](subData[1])
        ϵFunction = x -> convert_n2ϵ(nkFunc(x))
        permittivity = FunctionPermittivity(ϵFunction)
    end

    wavenumberRange = (WavenumberByλ₀(minλ₀), WavenumberByλ₀(maxλ₀))

    return Material(permittivity; wavenumberRange = wavenumberRange )
end

# Get the minimum and maximum wavelengths (in terms of microns) allowed by the model
function getWavelengthEnds(subData::Dict{Any,Any})
    wavelengthEnds = split(subData["wavelength_range"], " ")
    minλ₀ = parse(Float64, wavelengthEnds[1])*μm
    maxλ₀ = parse(Float64, wavelengthEnds[2])*μm
    return minλ₀, maxλ₀
end

# Parse the coefficients string into a float vector
function getCoefficients(subData::Dict{Any,Any})
    coefficientStr = split(subData["coefficients"], " ")
    return [parse(Float64, coefficient) for coefficient in coefficientStr]
end

function sortParallelVectors(xValues, yValues)
    dat = sort(collect(zip(xValues, yValues)); by=first)
    xValuesSorted = [dat[i][1] for i in 1:length(dat)]
    yValuesSorted = [dat[i][2] for i in 1:length(dat)]
    return xValuesSorted, yValuesSorted
end

# Create function giving complex n for a given wavenumber
function parseTabulatednk(subData::Dict{Any,Any})
    rawData = subData["data"]
    lines = split(rawData, "\n")
    nLines = length(lines)-1

    # Create parallel vectors defining λ₀, n, k
    λ₀Values = Vector{Float64}(undef,nLines)
    nRealValues = Vector{Float64}(undef,nLines)
    nImagValues = Vector{Float64}(undef,nLines)
    for iLine in 1:nLines
        values = split(lines[iLine], " ")
        λ₀Values[iLine] = parse(Float64,values[1])
        nRealValues[iLine] = parse(Float64,values[2])
        nImagValues[iLine] = parse(Float64,values[3])
    end

    # Convert to complex n, and put in SI units.
    nValues = nRealValues .+ nImagValues*1im
    λ₀Values = λ₀Values * μm

    # Sort
    λ₀ValuesSorted, nValuesSorted = sortParallelVectors(λ₀Values, nValues)


    # Create function to interpolate at a given wavenumber
    λnInterpolation = LinearInterpolation(λ₀ValuesSorted, nValuesSorted)
    function interpolateFunction(wavenumber::Wavenumber)::ComplexF64
        return λnInterpolation(getλ₀(wavenumber))
    end

    minλ₀ = minimum(λ₀ValuesSorted)
    maxλ₀ = maximum(λ₀ValuesSorted)

    return interpolateFunction, minλ₀, maxλ₀
end

# Create a function giving real n or real k for a given wavenumber
function parseSingularTabulatedx(subData::Dict{Any,Any})  # x can refer to n or k
    rawData = subData["data"]
    lines = split(rawData, "\n")
    nLines = length(lines)-1

    # Create vectors defining the wavelength and (n or k)
    λ₀Values = Vector{Float64}(undef,nLines)
    xValues = Vector{Float64}(undef,nLines)
    for iLine in 1:nLines
        values = split(lines[iLine], " ")
        λ₀Values[iLine] = parse(Float64,values[1])
        xValues[iLine] = parse(Float64,values[2])
    end

    # Convert to SI units
    λ₀Values = λ₀Values * μm

    # Sort
    λ₀ValuesSorted, xValuesSorted = sortParallelVectors(λ₀Values, xValues)

    # Create function for getting (n or k) at a given wavenumber
    λxInterpolation = LinearInterpolation(λ₀ValuesSorted, xValuesSorted) # create interpolation function
    function interpolateFunction(wavenumber::Wavenumber)::ComplexF64
        return λxInterpolation(getλ₀(wavenumber))
    end

    minλ₀ = minimum(λ₀ValuesSorted)
    maxλ₀ = maximum(λ₀ValuesSorted)

    return interpolateFunction, minλ₀, maxλ₀
end

# Formula1 = Sellmeier: n^2-1 = C1+ C2*λ²/(λ²-C3²) + C4*λ²/(λ²-C5²) + C4*λ²/(λ²-C5²) +...
function parseFormula1(subData::Dict{Any,Any})

    minλ₀, maxλ₀ = getWavelengthEnds(subData)
    coefficients = getCoefficients(subData)

    nPairs = half(length(coefficients)-1)

    function formula1(wavenumber::Wavenumber)::Float64
        local λ₀μm = getλ₀(wavenumber)/μm
        local n² = 1 + coefficients[1]
        local offset = 1
        for iPairs in 1:nPairs
            local C = coefficients[(offset+1):(offset+2)]
            n² += C[1]*λ₀μm^2 / (λ₀μm^2 - C[2]^2)
            offset += 2
        end
        return sqrt(n²)
    end

    return formula1, minλ₀, maxλ₀
end

# Formula2 = Sellmeier-2: n^2-1 = C1+ C2*λ²/(λ²-C3) + C4*λ²/(λ²-C5) + C4*λ²/(λ²-C5) +...
function parseFormula2(subData::Dict{Any,Any})

    minλ₀, maxλ₀ = getWavelengthEnds(subData)
    coefficients = getCoefficients(subData)

    nPairs = half(length(coefficients)-1)

    function formula2(wavenumber::Wavenumber)::Float64
        local λ₀μm = getλ₀(wavenumber)/μm
        local n² = 1 + coefficients[1]
        local offset = 1
        for iPairs in 1:nPairs
            local C = coefficients[(offset+1):(offset+2)]
            n² += C[1]*λ₀μm^2 / (λ₀μm^2 - C[2])
            offset += 2
        end
        return sqrt(n²)
    end

    return formula2, minλ₀, maxλ₀
end

# Formula 3 = Polynomial: n^2=C_1+C_2 λ^(C_3 )+C_4 λ^(C_5 )+C_6 λ^(C_7 )+C_8 λ^(C_9 )+C_10 λ^(C_11 )+C_12 λ^(C_13 )+C_14 λ^(C_15 )+C_16 λ^(C_17 )
function parseFormula3(subData::Dict{Any,Any})

    minλ₀, maxλ₀ = getWavelengthEnds(subData)
    coefficients = getCoefficients(subData)

    nPairs = half(length(coefficients)-1)

    function formula3(wavenumber::Wavenumber)::Float64
        local λ₀μm = getλ₀(wavenumber)/μm
        local n² = coefficients[1]
        local offset = 1
        for iPairs in 1:nPairs
            local C = coefficients[(offset+1):(offset+2)]
            n² += C[1]*λ₀μm^C[2]
            offset += 2
        end
        return sqrt(n²)
    end

    return formula3, minλ₀, maxλ₀
end

# Formula 4 = Refractiveindex.INFO type: n²=C_1+(C_2 λ^(C_3 ))/(λ^2-〖C_4〗^(C_5 ) )+(C_6 λ^(C_7 ))/(λ^2-〖C_8〗^(C_9 ) )+C_10 λ^(C_11 )+C_12 λ^(C_13 )+C_14 λ^(C_15 )+C_16 λ^(C_17 ) +...
function parseFormula4(subData::Dict{Any,Any})

    minλ₀, maxλ₀ = getWavelengthEnds(subData)
    coefficients = getCoefficients(subData)

    nQuads = 2 # There are always 2 sets of quads.  Variable number of pairs.
    nPairs = half(length(coefficients)-nQuads*4-1)

    function formula4(wavenumber::Wavenumber)::Float64
        local λ₀μm = getλ₀(wavenumber)/μm
        local n² = coefficients[1]
        local offset = 1
        for iQuads in 1:nQuads
            C = coefficients[(offset+1):(offset+4)]
            n² += (C[1] * λ₀μm^C[2]) / (λ₀μm^2 - C[3]^C[4])
            offset += 4
        end
        for iPairs in 1:nPairs
            C = coefficients[(offset+1):(offset+2)]
            n² += C[1] * λ₀μm^C[2]
            offset += 2
        end
        return sqrt(n²)
    end

    return formula4, minλ₀, maxλ₀
end

# Formula5 = Cauchy: =C_1+C_2 λ^(C_3 )+C_4 λ^(C_5 )+C_6 λ^(C_7 )+C_8 λ^(C_9 )+...
function parseFormula5(subData::Dict{Any,Any})

    minλ₀, maxλ₀ = getWavelengthEnds(subData)
    coefficients = getCoefficients(subData)

    nPairs = half(length(coefficients)-1)

    function formula5(wavenumber::Wavenumber)::Float64
        local λ₀μm = getλ₀(wavenumber)/μm
        local n = coefficients[1]
        local offset = 1
        for iPairs in 1:nPairs
            C = coefficients[(offset+1):(offset+2)]
            n += C[1]*λ₀μm^C[2]
            offset += 2
        end
        return n
    end

    return formula5, minλ₀, maxλ₀
end

# Formula6 = Gases: n-1 = C1 + C2/(C3-λ⁻²) + C4/(C5-λ⁻²) + C6/(C7-λ⁻²) +...
function parseFormula6(subData::Dict{Any,Any})

    minλ₀, maxλ₀ = getWavelengthEnds(subData)
    coefficients = getCoefficients(subData)

    nPairs = half(length(coefficients)-1)

    function formula6(wavenumber::Wavenumber)::Float64
        local λ₀μm = getλ₀(wavenumber)/μm
        local n = 1 + coefficients[1]
        local offset = 1
        for iPairs in 1:nPairs
            C = coefficients[(offset+1):(offset+2)]
            n += C[1] / (C[2] - λ₀μm^(-2))
            offset += 2
        end
        return n
    end

    return formula6, minλ₀, maxλ₀
end

# Formula 7 = Herzberger: n=C_1+C_2/(λ^2-0.028)+C_3 (1/(λ^2-0.028))^2+C_4 λ^2+C_5 λ^4+C_6 λ^6 ...
function parseFormula7(subData::Dict{Any,Any})

    minλ₀, maxλ₀ = getWavelengthEnds(subData)
    coefficients = getCoefficients(subData)

    # There are 3 special coefficients.  The rest are independent
    nSingles = length(coefficients)-3

    function formula7(wavenumber::Wavenumber)::Float64
        local λ₀μm = getλ₀(wavenumber)/μm
        local n = coefficients[1] + coefficients[2]/(λ₀μm^2-0.028) + coefficients[3]*(1/(λ₀μm^2-0.028))^2
        local offset = 3
        for iSingles in 1:nSingles
            C = coefficients[(offset+1)]
            n += C * λ₀μm^(2*iSingles)
            offset += 1
        end
        return n
    end

    return formula7, minλ₀, maxλ₀
end




# Formula 8 = Retro: (n^2-1)/(n^2+2)=C_1+(C_2 λ^2)/(λ^2-C_3 )+C_4 λ^2
function parseFormula8(subData::Dict{Any,Any})

    minλ₀, maxλ₀ = getWavelengthEnds(subData)
    coefficients = getCoefficients(subData)

    function formula8(wavenumber::Wavenumber)::Float64
        local λ₀μm = getλ₀(wavenumber)/μm
        local C = coefficients

        # n²m1Overn²p2 := (n^2-1)/(n^2+2)
        local n²m1Overn²p2 = C[1] + (C[2]*λ₀μm^2)/(λ₀μm^2-C[3]) + (C[4]*λ₀μm^2)

        local n² = (n²m1Overn²p2*2 + 1) / (1 - n²m1Overn²p2)

        return sqrt(n²)
    end

    return formula8, minλ₀, maxλ₀
end


# Formula 9 = Exotic: n^2=C_1+C_2/(λ^2-C_3 )+(C_4 (λ-C_5))/(〖(λ-C_5)〗^2+C_6 )
function parseFormula9(subData::Dict{Any,Any})

    minλ₀, maxλ₀ = getWavelengthEnds(subData)
    coefficients = getCoefficients(subData)

    function formula9(wavenumber::Wavenumber)::Float64
        local λ₀μm = getλ₀(wavenumber)/μm
        local C = coefficients

        local n² = C[1] + C[2]/(λ₀μm^2-C[3]) + (C[4]*(λ₀μm-C[5]))/( (λ₀μm-C[5])^2 + C[6])

        return sqrt(n²)
    end

    return formula9, minλ₀, maxλ₀
end


# Dictionary of functions, so the type-string can be matched to the analysis function
const global MATERIALDATABASEFUNCTIONS = Dict{String,Function}(
    "tabulated nk" => parseTabulatednk,
    "tabulated n" => parseSingularTabulatedx, # k is assumed to be 0
    "tabulated k" => parseSingularTabulatedx,
    "formula 1" => parseFormula1,
    "formula 2" => parseFormula2,
    "formula 3" => parseFormula3,
    "formula 4" => parseFormula4,
    "formula 5" => parseFormula5,
    "formula 6" => parseFormula6,
    "formula 7" => parseFormula7,
    "formula 8" => parseFormula8,
    "formula 9" => parseFormula9
    )
