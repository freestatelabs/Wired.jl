""" Defines output fields (3D spatial positions) for Wired.jl
    (c) 2024 ryan@freestatelabs

    Includes:
        - line (start, stop, number of points)
        - plane (corner 1 point, corner 2 point, number of points in both directions)
        - mesh
"""

"""
    abstract type Field

Represents a collection of locations in 3D space at which output data will be requested
"""
abstract type Field end 


"""
    line(start, stop, N)

Represents points along a line in 3D space

# Fields 
- `start::Vector{Float64}`: 3-length start point of the line 
- `stop::Vector{Float64}`: 3-length end point of the line 
- `N::Integer`: number of points along the line 
- `nodes::Matrix{Float64}: Nx3 `Matrix` of locations for calculations
"""
struct Line <: Field 

    start::Vector{Float64}
    stop::Vector{Float64}
    N::Integer
    nodes::Matrix{Float64}


    function Line(start::Vector{<:Number}, stop::Vector{<:Number}, N::Integer)

        nodes = zeros(N,3)
        nodes[1,:] = start 
        nodes[end,:] = stop 

        u = (stop .- start)/ (N - 1.0)

        for i in range(2, size(nodes)[1] - 1)
            nodes[i,:] = nodes[i-1,:] .+ u 
        end 

        new(float.(start), float.(stop), N, nodes)
    end
end 

    