""" Custom type definition for Wired.jl
    (c) 2024 ryan@freestatelabs
"""

abstract type Source end 


"""
    struct Wire <: Source 

Represents a finite wire segment in 3D space 

The wire segment has a finite radius, such that numerical singularities are 
prevented: the magnetic field generated by the wire segment has a finite (maximum)
value on the surface at r=R, and decays linearly to zero at the center with r=0.

# Fields 
- `a0::Vector{Float64}`: XYZ coordinates of the start of the wire vector
- `a1::Vector{Float64}`: XYZ coordinates of the end of the wire vector 
- `I::Float64`: total current in the filament
- `R::Float64`: radius of the wire segment
"""
struct Wire <: Source
    a0::SVector{3, Float32}
    a1::SVector{3, Float32}
    I::Float32
    R::Float32

    # function Wire(a0::SVector{<:Real}, a1::SVector{<:Real}, I::Real, R::Real)
    #     new(float.(a0), float.(a1), float(I), float(R))
    # end
end


abstract type Ring <: Source end 
"""
    struct CircularRing <: Ring

A circular current-carrying solid conducting ring with circular cross-section 

The ring has a defined cross-section, such that numerical singularities are prevented 
and the magnetic field can be accurately calculated within the ring itself.

# Fields 
- `name::String`: describes the filament
- `H::Float64`: distance along the Z-axis of the ring from the coordinate 
    system centroid to the ring centroid
- `R::Float64`: major radius of the ring
- `r::Float64`: minor radius of the ring (cross-section radius)
- `I::Float64`: current in the conducting ring; sign convention follows right-hand 
    rule about primary axis.
"""
struct CircularRing <: Ring 

    name::String
    H::Float64
    R::Float64 
    r::Float64 
    I::Float64 

    function CircularRing(name::AbstractString, H::Number, R::Number, r::Number, I::Number)

        new(name, float(H), float(R), float(r), float(I))
    end
end


"""
    struct RectangularRing <: Ring  

A circular current-carrying solid conducting ring with rectangular cross-section 

The ring has a defined cross-section, such that numerical singularities are prevented 
and the magnetic field can be accurately calculated within the ring itself.

# Fields 
- `name::String`: describes the filament
- `H::Float64`: distance along the Z-axis of the ring from the coordinate 
    system centroid to the ring centroid
- `R::Float64`: major radius of the ring
- `w::Float64`: width of the ring cross-section 
- `h::Float64`: height of the ring-cross-section
- `I::Float64`: current in the conducting ring; sign convention follows right-hand 
    rule about primary axis.
"""
struct RectangularRing <: Ring 

    name::AbstractString
    H::Number 
    R::Number 
    w::Number 
    h::Number
    I::Number 

    function RectangularRing(name::AbstractString, H::Number, R::Number, w::Number, h::Number, I::Number)

        new(name, float(H), float(R), float(w), float(h), float(I))
    end
end
