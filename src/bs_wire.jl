""" Biot-Savart law integration for finite wire segments 
    (c) 2024 ryan@freestatelabs 

    Defines 2 functions with 6 methods:
    - biotsavart!(): for actually calculating the Bfield 
    - 
"""


"""
    biotsavart!(nodes::AbstractArray{Float64}, wires::Vector{Wire})

Calculate the magnetic flux density generated by a series of current-carrying 
wire segments. 

Modifies an existing output array for the B-field in-place. Performs a current 
density correction for node points within the radius of the wire segment. Zeroes 
out singularity points.

"""
@views function biotsavart!(B::AbstractArray, nodes::AbstractArray{Float32}, wires::AbstractArray{Wire})

    # Prefer initializing to zero rather than undef because there's some sort of 
    # assignment stability issue?
    Nn = size(nodes)[1]
    b = zeros(Nn,3); c = zeros(Nn,3); cxa = zeros(Nn,3); 
    rp = zeros(Nn); rm = zeros(Nn)
    a = zeros(3)
    norm_cxa = zeros(Nn); dot_ac = zeros(Nn); 
    norm_c = zeros(Nn); dot_ab = zeros(Nn) 
    norm_b = zeros(Nn) 
    d = 0.0
    e = zeros(Nn)

    # Calculate the effect of each source on all nodes 
    # Linearly superimpose (sum) that effect from all sources
    for wire in wires

        d = mu0 * wire.I / (4pi)
        a .= wire.a1 .- wire.a0
        b .= wire.a0' .- nodes 
        c .= wire.a1' .- nodes
        R = wire.R

        crossrows!(cxa, c, a)
        normrows!(norm_cxa, cxa)
        dotrows!(dot_ac, c, a) 
        dotrows!(dot_ab, b, a) 
        normrows!(norm_c, c) 
        normrows!(norm_b, b)
        rp .= norm_cxa ./ norm(a)

        # Saving this as its own vector reduces allocations significantly, but 
        # doesn't have an effect on execution time
        e .= d .* (norm_cxa.^(-2)) .* (dot_ac./norm_c .- dot_ab./norm_b)
        multrows!(cxa, e) 

        # Reduce the current density if inside the conductor radius
        map!(x -> x < R ? (x^2 / R^2) : 1.0, rm, rp)

        B .+= cxa .* rm
        map!(x -> isnan(x) ? 0.0 : x, B, B)       # adds 0.4 ms
    end

    return B
end


"""
    biotsavart(nodes::AbstractArray{Float64}, wires::Vector{Wire})

Calculate the magnetic flux density at points in 3D space generated by a series 
of finite wire segments
"""
function biotsavart(nodes::AbstractArray{Float32}, wires::AbstractArray{Wire})
    
    B = zeros(size(nodes))
    biotsavart!(B, nodes, wires)

    return B 
end 


"""
    bfield(nodes::AbstractArray, wires::Vector{Wire}; Nt::Integer=0)


"""
function bfield(nodes::AbstractArray, wires::Vector{Wire}; Nt::Integer=0)

    Ns = length(wires)
    if Nt == 0 
        # Default is to use all available threads
        Nt = Threads.nthreads()
    end

    # Spawn a new task for each thread by splitting up the source array
    tasks = Vector{Task}(undef, Nt)
    for it = 1:Nt 
        @views tasks[it] = Threads.@spawn biotsavart(nodes, wires[threadindices(it, Nt, Ns)])
    end
    
    # Get the result from each calculation and add it to the output array 
    B = zeros(Float32, size(nodes))
    for it = 1:Nt 
        B .+= fetch(tasks[it]) 
    end 

    return B
end 