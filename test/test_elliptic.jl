""" Test and benchmark custom elliptic integral functions
"""

using Wired
using Elliptic, Printf, BenchmarkTools

# Wired.errmax = 1e-16
verbose = false
benchmark = false
errtol = 1e-6


function test_elliptic()

    println("Testing elliptic functions")
    N = 1000 
    k2 = [x for x in range(0, 0.999, N)]

    errE = abs.(ellipE.(k2; errmax=Wired.errmax) .- Elliptic.E.(k2))
    errK = abs.(ellipK.(k2; errmax=Wired.errmax) .- Elliptic.K.(k2))

    if verbose
        println("Max relative error: ")
        @printf "E() functions: %.3e\n" maximum(errE)
        @printf "K() functions: %.3e\n" maximum(errK)
    
        println("Avg relative error: ")
        @printf "E() functions: %.3e\n" mean(errE)
        @printf "K() functions: %.3e\n" mean(errK)
    end

    if (maximum(errE) <= errtol) && (maximum(errK) <= errtol)
        return true 
    else
        return false
    end
end


function benchmark_ellipticWired(k2) 
    ellipE.(k2; Wired.errmax);
    ellipK.(k2; Wired.errmax);
end


function benchmark_Elliptic(k2)
    Elliptic.E.(k2);
    Elliptic.K.(k2);
end


if benchmark
    println("Benchmarking Elliptic.jl K() and E() functions: ")
    @btime benchmark_Elliptic($k2)

    println("Benchmarking Wired ellipE and ellipK functions: ")
    @btime benchmark_ellipticWired($k2)
end



