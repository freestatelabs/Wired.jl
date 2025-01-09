using Wired
using Test, BenchmarkTools

@testset "Wired.jl" begin

    include("test_elliptic.jl")
    @test test_elliptic()

    include("test_wire.jl")
    println("Using Julia kernel")
    Wired.kernel = "julia"
    @test testwire1()
    @test testwire2()
    @test testwire3()
    println("Using C kernel")
    Wired.kernel = "c"
    Wired.precision = Float32
    @test testwire1()
    @test testwire2()
    @test testwire3()
    Wired.precision = Float64

    include("test_rings.jl")
    @test testring_circular()
    @test testring_rectangular()

    include("validation.jl")
    @test test_validation()

    include("test_fields.jl")
    @test test_line()

end
