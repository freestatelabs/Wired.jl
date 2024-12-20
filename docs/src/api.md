# API

Documentation for `Wired.jl`'s public-facing API. 

```@meta
CurrentModule = Wired 
```

## Types
Define custom types.
```@docs
Mesh
Wire
CircularRing
RectangularRing
```

## Mathematical Functions
Define high-performance version of elliptic integral functions optimized for use
with this code.

```@docs
ellipK
ellipE
ellipKE
```

## Processing Routines

```@docs
makewires
makecircrings
``` 

## File I/O
Define functions for loading and save data to/from files.

```@docs
loadmesh 
loadrings
```

## Output Fields 
Define functions that create basic sets of points at which to calculate fields.

```@docs 
Line
```

## Biot-Savart Solvers
Define functions that calculate the magnetic flux density due to defined 
current sources. 

```@docs
biotsavart!
```
