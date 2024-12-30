""" Wired.jl 
    (c) 2024 ryan@freestatelabs

A high-performance Biot-Savart law integrator, in Julia.

This code is released under an open-source MIT license. You may use this code 
for any lawful purpose. However, attribution is kindly requested. 

All units are MKS: length in METER, current in AMP, current density in AMP/METER^2, 
    magnetic flux density in TESLA.
"""

module Wired

const version = 1.0

using LinearAlgebra: norm, dot, det
using DelimitedFiles
using Logging 
import Elliptic
using StaticArrays

# Permeability of free space
const mu0 = 4pi * (1e-7)
export mu0

# Define floating-point precision and error tolerance of elliptic functions
precision = Float64
errmax = 1e-8

include("sources.jl")
export Source, Wire, Ring, CircularRing, RectangularRing

include("fields.jl")
export Line

include("mesh.jl")
export Mesh

include("utils.jl")
export ellipK, ellipE, ellipKE, crossrows!, normrows!, multrows!, dotrows!

include("processing.jl")
export makewires, makecircrings

include("io.jl")
export loadmesh, savemesh, loadrings, saverings, loadwires, savewires

include("bs_ring.jl")
include("bs_wire.jl")
export bfield

include("solve.jl")

include("lorentz.jl")
export lorentz, netload

end # module
