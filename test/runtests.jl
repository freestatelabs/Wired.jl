using Wired
using Test, BenchmarkTools

@testset "Wired.jl" begin

    include("test_elliptic.jl")
    @test test_elliptic()

    include("test_wire.jl")
    include("test_rings.jl")
    println("SETTING PRECISION TO DOUBLE")
    Wired.precision = Float64
    println("USING JULIA KERNEL")
    Wired.kernel = "julia"
    @test testwire1()
    @test testwire2()
    @test testwire3()
    @test testring_circular()
    @test testring_rectangular()
    println("USING C KERNEL")
    Wired.kernel = "c"
    @test testwire1()
    @test testwire2()
    @test testwire3()
    # @test testring_circular()
    # @test testring_rectangular()
    println("SETTING PRECISION TO SINGLE")
    Wired.precision = Float32
    println("USING JULIA KERNEL")
    Wired.kernel = "julia"
    @test testwire1()
    @test testwire2()
    @test testwire3()
    @test testring_circular()
    @test testring_rectangular()
    println("USING C KERNEL")
    Wired.kernel = "c"
    @test testwire1()
    @test testwire2()
    @test testwire3()
    # @test testring_circular()
    # @test testring_rectangular()
    Wired.precision = Float64


    include("validation.jl")
    # @test test_validation()

    include("test_fields.jl")
    @test test_line()

end
