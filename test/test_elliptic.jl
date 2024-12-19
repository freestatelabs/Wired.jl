using Wired

using Elliptic, Printf

N = 1000 
k2 = [x for x in range(0, 0.999, N)]
errE = zeros(N)
errK = zeros(N)

errE .= abs.(ellipE.(k2) .- Elliptic.E.(k2))
errK .= abs.(ellipK.(k2) .- Elliptic.K.(k2))


function testellipticWired(k2) 
    
    ellipE.(k2);
    ellipK.(k2);
end

function testElliptic(k2)

    Elliptic.E.(k2);
    Elliptic.K.(k2);

end

println("Testing Elliptic.jl K() and E() functions: ")
@btime testElliptic($k2)

println("Testing Wired ellipE and ellipK functions: ")
@btime testellipticWired($k2)

println("Max relative error: ")
@printf "E() functions: %.3e\n" maximum(abs.(errE))
@printf "K() functions: %.3e\n" maximum(abs.(errK))

