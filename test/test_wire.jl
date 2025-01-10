using Wired


function testwire1()
    # Test calculation of Bfield generated by infinite Wire with infinitesimal cross-section

    println("Testing Wire - Infinitesimal Cross-Section")

    Iwire = 1000
    r = 1.0
    wire = Wire([0,0,-10000],[0,0,10000],Iwire,0)
    loc = [r 0 0]

    B = bfield(loc, [wire])

    if isapprox(B[2], mu0*Iwire/(2*pi*r),atol=1e-8)
        return true 
    else 
        return false 
    end 

end 


function testwire2()
    # Test calculation of Bfield generated by infinite Wire with finite cross-section

    println("Testing Wire - Finite Cross-Section")

    Iwire = 1000
    r = 1.0
    R = 1.0
    wire = Wire([0,0,-10000],[0,0,10000],Iwire,R)
    loc = [r 0 0]

    B = bfield(loc, [wire])

    if isapprox(B[2], mu0*Iwire*r/(2*pi*R^2),atol=1e-8)
        return true 
    else 
        return false 
    end 

end 

function analytical_wire(I, R, r)
    B = zeros(length(r))

    for i in eachindex(r)
        if r[i] < R 
            B[i] = mu0*I*r[i]/(2*pi*R^2) 
        else
            B[i] = mu0*I/(2*pi*r[i])
        end
    end

    return B 
end

function testwire3()
    println("Testing Wire - Finite Cross-Section, N=100")

    Iwire = 1000
    r = 5.0
    R = 1.0
    wire = Wire([0,0,-10000],[0,0,10000],Iwire,R)
    nodes = Line([0,0,0],[r,0,0],100).nodes

    By = bfield(nodes, [wire])[:,2]
    By_ = analytical_wire(Iwire, R, nodes[:,1])

    delta = By .- By_ 
    maxabserr = maximum(abs.(delta))

    if maxabserr < 1e-8
        return true 
    else 
        return false 
    end 
end

function plotdiff()
    Wired.precision=Float32
    Iwire = 1000
    r = 5.0
    R = 1.0
    wire = Wire([0,0,-10000],[0,0,10000],Iwire,R)
    nodes = Line([0,0,0],[r,0,0],100).nodes
    x = nodes[:,1]


    Wired.kernel="julia"
    By = bfield(nodes, [wire])[:,2]
    Wired.kernel="c"
    Byc = bfield(nodes, [wire])[:,2]
    By_ = analytical_wire(Iwire, R, nodes[:,1])
    plot(x, By, label="julia")
    plot!(x,Byc,label="c")
    plot!(x,By_,label="analytical")
end
# plotdiff()