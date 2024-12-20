""" Build documentation for Wired.jl
"""

using Documenter
# using LiveServer
using Wired

pages = [
    "Home" => "index.md", 
    "Basic Usage" => "usage.md",
    "Theory" => "theory.md", 
    "Validation" => "validation.md", 
    "Benchmarking" => "benchmarking.md", 
    "API" => "api.md"
    ]


makedocs(
    sitename="Wired.jl", 
    pages = pages)

# serve(dir="Wired/docs/build")

deploydocs(
    repo = "github.com/freestatelabs/Wired.jl.git",
)