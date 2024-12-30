""" Wired.jl
    Defines output fields (3D spatial positions)
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
struct Line{T<:Real} <: Field 

    start::Vector{T}
    stop::Vector{T}
    N::Integer
    nodes::Matrix{T}

    function Line{T}(start::Vector{<:Real}, stop::Vector{<:Real}, N::Integer) where T<:Real

        nodes = zeros(T, N, 3)
        nodes[1,:] = start 
        nodes[end,:] = stop 

        u = (stop .- start)/ (N - 1.0)

        for i in range(2, size(nodes)[1] - 1)
            nodes[i,:] = nodes[i-1,:] .+ u 
        end 

        new(convert.(T, start), convert.(T, stop), N, convert.(T, nodes))
    end

    # Convenience constructor for using Wired.precision
    function Line(start::Vector{<:Real}, stop::Vector{<:Real}, N::Integer)
        Line{precision}(start, stop, N)
    end
end