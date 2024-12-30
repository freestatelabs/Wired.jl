""" Performance benchmarks for Wired.jl
    (c) 2024 ryan@freestatelabs
"""

using Wired
using BenchmarkTools, Printf

run_wirebenchmarks = true
run_ringbenchmarks = false
make_plots = false
Wired.precision = Float64
Nvals = [1000, 2000, 4000, 8000, 10000]
Nts = [1, 2, 4, 6]

if make_plots	
	# Don't want to make Plots a dependency
	using Plots 
end

function createwireproblem(N)
    # Create a Wire source problem that results in N^2 Biot-Savart computations
    # Nnodes = Nwires = N 

    nodes = rand(Wired.precision, N, 3) 
    wires = Vector{Wired.Wire{Wired.precision}}(undef, N) 

    for i in range(1, N) 
        wires[i] = Wire(rand(3), rand(3), randn(), randn())
    end 

    return nodes, wires
end 


function createringproblem(N)
    # Create a Wire source problem that results in N^2 Biot-Savart computations
    # Nnodes = Nrings = N 

    nodes = rand(Wired.precision, N, 3) 
    rings = Vector{Wired.CircularRing{Wired.precision}}(undef, N) 

    for i in range(1, N) 
        rings[i] = CircularRing("", randn(), randn(), randn(), randn())
    end 

    return nodes, rings
end 


function runwirebenchmarks(Nvals::AbstractVector, Nts::AbstractVector)
    # Run Wire benchmarks for a series of N and Nt values

    # Each row corresponds to a value of N, each column is a value of Nt
    times = zeros(length(Nvals), length(Nts))
    println("Datatype set to: " *string(Wired.precision))  

    for i in range(1, length(Nvals))

        for j in range(1, length(Nts))
            nodes, wires = createwireproblem(Nvals[i])

            @printf "Running Wire benchmark with N=%i and Nt=%i\n" Nvals[i] Nts[j]
            trial = @benchmark bfield($nodes, $wires, Nt=$Nts[$j]);
            times[i,j] = median(trial).time
        end
    end

    return times ./ 1e6     # Time is in [ms]

end


function runringbenchmarks(Nvals::AbstractVector, Nts::AbstractVector)
    # Run Ring benchmarks for a series of N and Nt values

    # Each row corresponds to a value of N, each column is a value of Nt
    times = zeros(length(Nvals), length(Nts))
    println("Datatype set to: " *string(Wired.precision)) 

    for i in range(1, length(Nvals))

        for j in range(1, length(Nts))
            nodes, rings = createringproblem(Nvals[i])

            @printf "Running Ring benchmark with N=%i and Nt=%i\n" Nvals[i] Nts[j]
            trial = @benchmark bfield($nodes, $rings, Nt=$Nts[$j]);
            times[i,j] = median(trial).time
        end
    end

    return times ./ 1e6     # Time is in [ms]

end


function plotbenchmarks(Nvals, Nts, times, sourcename)
    p = plot()
    for i in range(1, length(Nts))
        plot!(p, Nvals, times[:,i], label="Nt="*string(Nts[i]))
    end
    xlabel!("N")
    ylabel!("Computation Time [ms]")
    title!(sourcename*" Source - Benchmarks\n(on 2021 M1 MacBook Pro)")

    return p

end


if run_wirebenchmarks

    times1 = runwirebenchmarks(Nvals, Nts)

    if make_plots
        p1 = plotbenchmarks(Nvals, Nts, times1, "Wire")
        display(p1)
    end
end 


if run_ringbenchmarks

    times2 = runringbenchmarks(Nvals, Nts)

    if make_plots
        p2 = plotbenchmarks(Nvals, Nts, times2, "Ring")
        display(p2)
    end
end


