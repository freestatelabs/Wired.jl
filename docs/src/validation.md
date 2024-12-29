# Validation

Code validation for `Wired.jl`.

## Example Problem 

Given a circular curent-carrying ring with a circular cross-section and the following parameters:
- `H = 0.0 m`
- `R = 2.0 m`
- `r = 0.1 m` 
- `I = 10 kA`

Calculate the magnetic field at the following locations:
- Along the Z-axis from (0,0,0) to (0,0,1) [m]
- Along the X-axis from (0,0,0) to (0,0,4) [m]

With the following methods
- Analytical (valid only for Z-axis)
- Using a `Wired.CircularRing` primitive 
- Using a finite element mesh

## Code

```julia
H = 0.0 	# height above xy plane
R = 2.0 	# major radius
I = 10e3 	# 10kA current
r = 0.1 	# minor radius (cross-section)

xaxis = Line([0,0,0], [4.0,0,0], 1000)
zaxis = Line([0,0,0], [0,0,1.0], 1000)  

analytical(z) = mu0 * I * (R^2) ./ (2 .* (R^2 .+ z.^2).^1.5)
ring = CircularRing("Ring", H, R, r, I) 

Bz_ring = bfield(zaxis.nodes, [ring])
Bz_analytical = analytical.(zaxis.nodes[:,3])
```

## Results

![Z-Axis Results](figs/validation-zaxis.svg)