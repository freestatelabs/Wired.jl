""" Finite element mesh operations for Wired.jl 
    (c) 2024 ryan@freestatelabs
"""

"""
    mutable struct Mesh  

Define a finite element mesh, which can act as a Source and a Field

# Fields 
- `nodes::AbstractArray`: Nx3 matrix of positions in 3D space (element centroids)
- `volumes::AbstractArray`: N-length vector of volume of each element 
- `Jdensity::AbstractArray`: Nx3 matrix of current density vectors at each element 
- `mu_r::Float64`: Relative magnetic permeability of the mesh
"""
mutable struct Mesh

    nodes::AbstractArray
    volumes::AbstractArray
    Jdensity::AbstractArray
    mu_r::Float64

    # Sometimes the current density is not defined
    function Mesh(nodes::AbstractArray, volumes::AbstractArray, mu_r::Number)
        
        new(nodes, volumes, [], float(mu_r))
    end 

    function Mesh(nodes::AbstractArray, volumes::AbstractArray, Jdensity::AbstractArray, mu_r::Number)

        new(nodes::AbstractArray, volumes::AbstractArray, Jdensity::AbstractArray, float(mu_r))
    end
end


"""
    tetvolume(cornernodes::AbstractArray)

Calculate the volume of a tetrahedral element using the spatial position of its 
four corner nodes. 

Ref: https://math.stackexchange.com/questions/3616760/how-to-calculate-the-volume-of-tetrahedron-given-by-4-points 
"""
function tetvolume(cornernodes::AbstractArray)


    return (1/6)*abs(det(hcat(cornernodes, ones(4))))
end