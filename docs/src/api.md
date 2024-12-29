# API

Documentation for `Wired.jl`'s public-facing API. 

```@meta
CurrentModule = Wired 
```

## Sources
Define custom source types.
```@docs
Source
Wire
Ring
CircularRing
RectangularRing
```

## Mathematical Functions
Define high-performance version of elliptic integral and linear algebra functions optimized for use
with this code.

```@docs
ellipK
ellipE
ellipKE
crossrows!
normrows!
multrows!
dotrows!
```

## Processing Routines

```@docs
makewires
makecircrings
threadindices
``` 

## File I/O
Define functions for loading and save data to/from files.

```@docs
loadmesh 
savemesh
loadrings
saverings
loadwires
savewires
```

## Output Fields 
Define functions that create basic sets of points at which to calculate fields.

```@docs 
Line
Mesh
```

## Biot-Savart Solvers
Define functions that calculate the magnetic flux density due to defined 
current sources. 

```@docs
bfield
```
