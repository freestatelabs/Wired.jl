using Wired
using Test

@testset "Wired.jl" begin

    include("test_wire.jl")
    @test testwire1()
    @test testwire2()

    include("test_rings.jl")
    @test testring_circular()
    @test testring_rectangular()

    include("validation.jl")
    @test test_validation()
end
