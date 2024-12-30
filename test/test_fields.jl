""" Wired.jl 
    Test field objects
"""

function test_line(start=[0,0,0], stop=[1,1,1], N=100)
    # Check that the line object instantiates correctly 
    # (currently incomplete)

    println("Testing Line")

    # Use Wired.precision
    line = Line(start, stop, N)

    if (size(line.nodes) == (N,3))
        return true 
    else
        return false 
    end
end