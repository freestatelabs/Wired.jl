""" Performance benchmarks for Wired.jl
    (c) 2024 ryan@freestatelabs
"""

using Wired
using Printf, StaticArrays, BenchmarkTools

# @printf "Running finite wire benchmarks...\n"

function createproblem(N)

    nodes = rand(Float32, N, 3) 
    wires = Vector{Wired.Wire}(undef, N) 

    for i in range(1, N) 
        a = @SVector [x for x in rand(Float32, 3)]
        b = @SVector [x for x in rand(Float32, 3)]

        wires[i] = Wire(a, b, rand(Float32, 1)[1], rand(Float32, 1)[1])
    end 

    return nodes, wires
end 


N = 1000 
Nt = 1
nodes, wires = createproblem(N) 

@btime bfield($nodes, $wires, Nt=Nt)