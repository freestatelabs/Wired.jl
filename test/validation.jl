""" Wired.jl
	Validation problem for the documentation; see
		freestatelabs.github.io/Wired/theory
"""

using Wired 
using Plots 


function test_validation(do_plot=false, savefigs=false)

	H = 0.0 
	R = 2.0 
	I = 10e3 
	r = 0.1

	xaxis = Line([0,0,0], [4.0,0,0], 1000)
	zaxis = Line([0,0,0], [0,0,1.0], 1000)  

	analytical(z) = mu0 * I * (R^2) ./ (2 .* (R^2 .+ z.^2).^1.5)

	ring = CircularRing("Ring", H, R, r, I) 
	# mesh = loadmesh("testdata/validation-mesh.csv")
	# meshwires = makewires(mesh)

	Bz_ring = bfield(zaxis.nodes, [ring])
	Bz_analytical = analytical.(zaxis.nodes[:,3])
	# Bz_mesh = bfield(zaxis.nodes, meshwires)

	By_ring = bfield(xaxis.nodes, [ring])
	By_analytical = analytical.(xaxis.nodes[:,1])

	if do_plot
		p = plot(zaxis.nodes[:,3], Bz_analytical, label="Analytical", color=:black)
		plot!(p, zaxis.nodes[1:50:1000,3], Bz_ring[1:50:1000,3], label="Ring Primitive", markershape=:circle, markercolor=:red, markeralpha=1.0, linetype=:scatter)
		xlabel!("Distance Along Z-axis [m]")
		ylabel!("Magnetic Flux Density (Bz) [T]")

		if savefigs
			savefig(p, "validation-zaxis.svg")
		end

		p2 = plot(xaxis.nodes[1:50:1000,1], By_ring[1:50:1000,1], label="Ring Primitive", markershape=:circle, markercolor=:red, markeralpha=1.0, linetype=:scatter)
		xlabel!("Distance Along X-axis [m]")
		ylabel!("Magnetic Flux Density (By) [T]")

		display(p2)
		if savefigs
			savefig(p2, "validation-xaxis.svg")
		end

	end

	if isapprox(Bz_analytical[1], Bz_ring[1,3], atol=1e-8)
		return true 
	else
		return false
	end 

end
