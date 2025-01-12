""" Build documentation for Wired.jl
"""

using Documenter
using Wired

pages = [
    "Home" => "index.md", 
    "Basic Usage" => "usage.md",
    "Theory" => "theory.md", 
    "Validation" => "validation.md", 
    "Benchmarking" => "benchmarking.md", 
    "C Kernel" => "kernel.md",
    "API" => "api.md"
    ]


makedocs(
    sitename="Wired.jl", 
    pages = pages)

# using LiveServer
# serve(dir="Wired/docs/build")

deploydocs(
    repo = "github.com/freestatelabs/Wired.jl.git",
)