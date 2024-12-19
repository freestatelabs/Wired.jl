using BenchmarkTools, Plots
using Wired

# # circ, rect = loadrings("testdata/testring.csv")
circ = CircularRing("name", 0, 1, 0.05, 1000)
rect = RectangularRing("name", 0, 1, 0.1, 0.1, 1000)

z = 0.
x0 = 0.
x1 = 2.
a = circ.R[1]
line = Line([x0,0,z], [x1,0,z], 1000) 

Bc = bfield(line.nodes, [circ])
Br = bfield(line.nodes, [rect]; Nmin=4)

# Baxis = length(circ.name)*(mu0 * circ.I[1] * a^2) / (2*(a^2 + z^2)^1.5)

p = plot(line.nodes[:,1], Bc[:,3], label="Circular")#, seriestype=:scatter)
plot!(line.nodes[:,1], Br[:,3], label="Rectangular")#, seriestype=:scatter)
display(p)

@info "Program completed successfully."

nodes = zeros(1000,3)
B = zeros(size(nodes))
circ = Vector{CircularRing}(undef,1000)
for i in range(1,length(circ)) 
    circ[i] = CircularRing("name", 0, 1, 0.05, 1000)
end

@benchmark bfield($nodes, $circ; errmax=1e-8, Nt=1)





