# Usage

Basic usage guide for `Wired.jl`.

*All dimensions are in MKS: length in meters, current in amperes, field in Tesla.*

## Fields 

All results are calculated on `Field` objects. At the most basic, a `Field` is a `Matrix` representing points in 3D space, i.e.:

```julia
julia> line = Line([0, 0, 0], [1, 0, 0], 11);

julia> line.nodes
11×3 Matrix{Float64}:
 0.0  0.0  0.0
 0.1  0.0  0.0
 0.2  0.0  0.0
 0.3  0.0  0.0
 0.4  0.0  0.0
 0.5  0.0  0.0
 0.6  0.0  0.0
 0.7  0.0  0.0
 0.8  0.0  0.0
 0.9  0.0  0.0
 1.0  0.0  0.0
```

## Sources 

### `Wire` Source
The most basic type of source is a `Wire`: a finite-length, circular conductor in 3D space. The `Wire` is completely defined by four parameters:
* The start location of the wire (3-length vector)
* The end location of the wire (3-length vector)
* The total current flowing in the wire
* The radius of the wire

```julia
julia> start = [0, 0, 0];
julia> stop = [3, 100, 0];
julia> current = 10;
julia> radius = 0.01;
julia> wire = Wire(start, stop, current, radius)
Wire([0.0, 0.0, 0.0], [3.0, 100.0, 0.0], 10.0, 0.01)
```


### `Ring` Sources
`Wired.jl` supports two kinds of `Ring` sources: `CircularRing`'s and `RectangularRing`'s. The major axis is the Z-axis. Each is defined by the following parameters:
- `H`: the height above the horizontal (XY) plane 
- `R`: the major radius, i.e. the radius from the central (Z) axis to the center of the ring cross-section
- `I`: the total current in the conductor 
- `r`: the minor radius (for a circular cross-section)
- `w`,`h`: the width and height (for a rectangular cross-section)

```julia
julia> H=0;    R=1.0;    I=1000;    r=0.1;

julia> ring = CircularRing("Circular Ring 1", H, R, r, I)
CircularRing("Circular Ring 1", 0.0, 1.0, 0.1, 1000.0)
```

```julia
julia> H=1.0;    R=2.0;    I=2000;    w=0.05;    h=0.10;

julia> ring = RectangularRing("Rectangular Ring 1", H, R, w, h, I)
RectangularRing("Rectangular Ring 1", 1.0, 2.0, 0.05, 0.1, 2000.0)
```


## Finite Element Meshes

`Wired.jl` provides operations for working with finite element meshes. All outputs are calculated at the centroid of the elements. The current density within the element is used to approximate a finite-length current-carrying `Wire` with circular cross-section (see [Finite Element Meshes]())

### Load a mesh from a file
`Wired.jl` assumes the file is comma-separated format, with the first three columns as x,y,z coordinates of the element centroid, the fourth column is the element volume, and the fifth through seventh columns are the elemental current density vectors.

```julia 
julia> mesh = loadmesh("testmesh.csv")

julia> wires = makewires(mesh)          # convert to Wire objects  

julia> B = bfield(mesh.nodes, wires)	# calculate self-field   
```