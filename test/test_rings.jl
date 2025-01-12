using Wired


function testring_circular()
    # Check that the magnetic field at the centroid of the ring is the same 
    # as the analytic expression

    println("Testing Ring - Circular")

    H = 0.0  
    R = 1.0
    r = 0.1 
    Iring = 1000
    circ = CircularRing("name", H, R, r, Iring)

    nodes = [0 0 0]
    B = bfield(nodes, [circ], errmax=1e-8)
    Bz = mu0 * Iring * R^2 / (2 * (R^2 + H^2)^1.5)

    if isapprox(B[3], Bz, rtol=1e-4)
        return true 
    else 
        return false 
    end
end


function testring_rectangular()
    # Check that the magnetic field at the centroid of the ring is the same 
    # as the analytic expression
    #
    # Note reduced precision in `rtol`: there's no analytic expression for this

    println("Testing Ring - Rectangular")

    H = 0.0  
    R = 10.0
    w = 0.1
    h = 0.1
    Iring = 1000

    rect = RectangularRing("name", H, R, w, h, 1000)
    nodes = [0.0 0.0 0.0] 
    B = bfield(nodes, [rect]; Nmin=10, errmax=1e-16)
    Bz = mu0 * Iring * R^2 / (2 * (R^2 + H^2)^1.5)

    if isapprox(B[3], Bz, rtol=1e-4)
        return true 
    else 
        return false 
    end

end


