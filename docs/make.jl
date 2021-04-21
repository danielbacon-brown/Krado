push!(LOAD_PATH,"../src/")
using Documenter, Krado

makedocs(sitename="Krado", modules=[Krado],
    pages = [
            "Index" => "index.md",
            "Lattice" => "Lattice.md",
        ]
    )

# deploydocs(
#     repo = "github.com/danielbacon-brown/Krado"
# )
