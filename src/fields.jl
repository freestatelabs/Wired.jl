""" Defines output fields (3D spatial positions) for Wired.jl
    (c) 2024 ryan@freestatelabs

    Includes:
        - line (start, stop, number of points)
        - plane (corner 1 point, corner 2 point, number of points in both directions)
        - mesh
"""

abstract type Field end 

"""
    line(start, stop, N)

Create points along a line in 3D space. 
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

    