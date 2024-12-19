""" Data processing routines for Wired.jl
    (c) 2024 ryan@freestatelabs
"""


"""
    makewires(mesh::Mesh)

Create finite Wire segment sources from a Mesh

# Arguments 
- `mesh::Mesh`: the finite element mesh to generate Wires from 
- `ratio=2.0`: the length/radius ratio for the wire segments. Set to 2.0 for perfect 
    cubes and 1.0 for tetrahedrons.

# Returns 
`Vector{Wire}` containing all Wire objects in the mesh
"""
function makewires(mesh::Mesh; ratio=2.0)

    # Assume each element is a perfect cube of side length L 
    # Convert this cube into a cylinder of length L and diameter L
    # This is also the half-step length for start/end points
    # Set L/r equal to `ratio`
    r = (mesh.volumes ./ (ratio*pi)).^(1.0/3.0)
    L = ratio.*r

    # Determine unit vectors from current densities
    Jmag = sqrt.(mesh.Jdensity[:,1].^2 + mesh.Jdensity[:,2].^2 + mesh.Jdensity[:,3].^2)
    u = mesh.Jdensity./Jmag
    map!(x -> isnan(x) ? 0.0 : x, u, u)     # Remove divide by zero issues
    
    # Determine magnitude of current flowing in this cylinder 
    I = pi .* (r.^2) .* Jmag 

    # Create start and end points 
    A = mesh.nodes .+ (L./2).*u 
    B = mesh.nodes .- (L./2).*u 

    wires = Vector{Wire}(undef, length(I))
    for i in range(1, length(I))
        wires[i] = Wire(A[i,:], B[i,:], I[i], r[i])
    end 

    return wires
end


"""
    function makecircrings(rect::RectangularRings; Nmin=2)

Filamentize Rings with rectangular cross-sections. 

For a cross-section with side lengths L1 and L2, with L1 < L2, the spacing `s` 
between filament centroids is given by s = L1 / Nmin. This spacing is used for 
creating filaments along both sides, such that N1 = Nmin and N2 ~= L2 / s (with 
integer rounding.)

# Arguments 
- `rect::Vector{RectangularRing}`: RectangularRing objects from which to create 
    CircularRing objects 
- 'Nmin=2`: Minimum number of CircularRings per side of the original cross-section. 

# Returns 
Vector{CircularRing}, which is of indeterminate length apriori.
"""
function makecircrings(rect::AbstractArray{RectangularRing}, Nmin=2)

    circ = Vector{CircularRing}(undef, 0)

    for i in eachindex(rect)

        # Determine number of filaments for each axis 
        if rect[i].w <= rect[i].h 
            Nw = Nmin 
            Nh = Nmin*floor(rect[i].h/rect[i].w)
            rfil = 0.5*rect[i].w / Nmin
        else 
            Nw = Nmin*floor(rect[i].w/rect[i].h)
            Nh = Nmin 
            rfil = 0.5*rect[i].h / Nmin
        end 

        Ifil = rect[i].I / (Nw * Nh)

        # Start from top left, work across and then down
        counter = 1
        # work in (u,v) coordinates, find the centroid of the first filament
        vc = rect[i].H + rect[i].h/2 - rfil 
        for j in range(1, Nh)

            uc = rect[i].R - rect[i].w/2 + rfil 
            for k in range(1, Nw)

                push!(circ, 
                    CircularRing(
                        string(rect[i].name, " ", string(counter)), 
                        vc, uc, rfil, Ifil
                    )
                )

                # Update horizontal position
                uc += 2 * rfil

                counter += 1
            end

            # Update vertical position 
            vc -= 2 * rfil

        end

    end

    return circ
end