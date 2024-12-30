using Wired
using Test, BenchmarkTools

@testset "Wired.jl" begin

    include("test_elliptic.jl")
    @test test_elliptic()

    include("test_wire.jl")
    @test testwire1()
    @test testwire2()

    include("test_rings.jl")
    @test testring_circular()
    @test testring_rectangular()

    include("validation.jl")
    @test test_validation()

    include("test_fields.jl")
    @test test_line()

end
