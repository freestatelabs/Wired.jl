""" Wired.jl 
    Biot-Savart law integration for finite wire segments 
"""

"""
    biotsavart!(B::AbstractArray{T}, nodes::AbstractArray{T}, 
                            wires::AbstractArray{Wire{T}}; mu_r=1.0) 

Calculate the magnetic flux density generated by a series of current-carrying 
wire segments. 

Modifies an existing output array for the B-field in-place. Performs a current 
density correction for node points within the radius of the wire segment. Zeroes 
out singularity points.
"""
@views function biotsavart!(B::AbstractArray{T}, nodes::AbstractArray{T}, 
                            wires::AbstractArray{Wire{T}}; mu_r=1.0) where T<:Real

    # Prefer initializing to zero rather than undef because there's some sort of 
    # assignment stability issue?
    Nn = size(nodes)[1]
    b = zeros(T, Nn,3); c = zeros(T, Nn,3); cxa = zeros(T, Nn, 3); 
    rp = zeros(T, Nn); rm = zeros(T, Nn)
    a = zeros(T, 3)
    norm_cxa = zeros(T, Nn); dot_ac = zeros(T, Nn); 
    norm_c = zeros(T, Nn); dot_ab = zeros(T, Nn) 
    norm_b = zeros(T, Nn) 
    d = 0.0
    e = zeros(T, Nn)

    # Calculate the effect of each source on all nodes 
    # Linearly superimpose (sum) that effect from all sources
    for wire in wires

        d = mu_r * mu0 * wire.I / (4pi)
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
function biotsavart(nodes::AbstractArray{T}, wires::AbstractArray{Wire{T}}; 
                    mu_r=1.0) where T<:Real
    
    B = zeros(T, size(nodes))
    biotsavart!(B, nodes, wires; mu_r=mu_r)

    return B 
end 


"""
    bfield(nodes::AbstractArray, wires::Vector{Wire}; 
            Nt::Integer=0, mu_r=1.0)

Calculate the B-field at a collection of points in 3D space, generated by a series of
finite-length `Wire` objects.

# Arguments
- `nodes::AbstractArray`: Nx3 `Matrix` containing (x,y,z) coordinates of points in 3D space
- `wires::Vector{Wire}`: `Wire` objects contributing to the magnetic field 
- `Nt::Integer`: number of threads to use for the calculation (default: all available threads)

# Returns
Nx3 `Matrix` containing magnetic flux density vectors at each of the points in 3D space represented by `nodes`
"""
function bfield(nodes::AbstractArray{T}, wires::Vector{Wire{S}}; 
                Nt::Integer=0, mu_r=1.0) where {T<:Real, S<:AbstractFloat}

    if T != S 
        nodes = convert.(S, nodes) 
    end
    Ns = length(wires)
    if Nt == 0 
        # Default is to use all available threads
        Nt = Threads.nthreads()
    elseif Nt > Threads.nthreads()
        println("Error. Number of threads specified is greater than available threads.")
    end

    # Spawn a new task for each thread by splitting up the source array
    tasks = Vector{Task}(undef, Nt)
    for it = 1:Nt 
        @views tasks[it] = Threads.@spawn biotsavart(nodes, wires[threadindices(it, Nt, Ns)]; mu_r=mu_r)
    end
    
    # Get the result from each calculation and add it to the output array 
    B = zeros(T, size(nodes))
    for it = 1:Nt 
        B .+= fetch(tasks[it]) 
    end 

    return B
end 