""" Wired.jl
	Validation problem for the documentation; see
		freestatelabs.github.io/Wired/theory
"""

using Wired 
using Printf
make_plots = false


if make_plots	
	# Don't want to make Plots a dependency
	using Plots 
end


function test_validation(make_plots=false, savefigs=false)

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

	if make_plots
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


function test_wireconvergence(Nwires; H=0, R=1, r=0.01, I=1e3)
	# Test how many wires are required to converge to ring solution 

	ring = CircularRing("Ring", H, R, r, I)

	# Discretize the ring into `Wire` segments
	wires = Vector{Wire{Float64}}(undef, Nwires)
	dtheta = (2pi)/(Nwires-1)
	theta = dtheta
	wires[1] = Wire{Float64}([R,0,0],[R*cos(theta), R*sin(theta), 0], I, r)
	for i in 2:Nwires 
		x0 = wires[i-1].a1[1] 
		y0 = wires[i-1].a1[2]
		wires[i] = Wire{Float64}([x0, y0, 0], [R*cos(theta),R*sin(theta),0], I, r)
		theta += dtheta
	end

	line = Line([0,0,0], [2,0,0], 1000)
	Bring = bfield(line.nodes, [ring], errmax=1e-16)
	Bwires = bfield(line.nodes, wires)

	@printf "Max/min Bring = %.3e/%.3e T\n" maximum(Bring) minimum(Bring)
	@printf "Max/min Bwires = %.3e/%.3e T\n" maximum(Bwires) minimum(Bwires)
	errmax = 100*(maximum(Bring) - maximum(Bwires)) / maximum(Bwires)
	errmin = 100*(minimum(Bring) - minimum(Bwires)) / minimum(Bwires)
	@printf "Rel error max/min = %.0f%%/%.0f%%" errmax errmin

	if make_plots
		p = plot()
		plot!(p, line.nodes[:,1], Bring[:,3], label="Bring")
		plot!(p, line.nodes[:,1], Bwires[:,3], label="Bwires")
		xlabel!("Distance from central axis [m]")
		ylabel!("Magnetic Flux Density (Bz) [T]")
		display(p)
	end

end
