using LatticeQCD
using Documenter

makedocs(
    sitename="LatticeQCD.jl",
    pages=[
        "Home" => "index.md",
        "Backend support" => "backends.md",
        "Migrating to v2" => "migration-v2.md",
    ],
)
