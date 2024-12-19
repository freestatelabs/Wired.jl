
using BenchmarkTools
using Wired

N = 1000
line = Line([0,0,0], [3,3,0], N)
wires = Vector{Wire}(undef, N) 
for i in eachindex(wires)
    wires[i] = Wire([0,0,-1e3],[0,0,1e3],5,1)
end
B = bfield(line.nodes, wires; Nt=0)
# plot(
#     scatter(
#         x = line.nodes[:,1], 
#         y = B[:,2]
#     )
# )

@benchmark bfield($line.nodes, $wires; Nt=0)

