""" Build documentation for Wired.jl
"""

using Documenter, LiveServer
cd("/Users/ryan/Github/Wired.jl/docs")
# push!(LOAD_PATH, ".../src/")
include("../src/Wired.jl")

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

serve(dir="build")