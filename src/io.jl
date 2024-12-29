""" I/O Functionality for Wired.jl
    (c) 2024 ryan@freestatelabs
"""


""" 
    loadmesh(fn::String)

Load a finite element mesh from a file

The first row is a header row. The subsequent rows are: 
Xcentroid, Ycentroid, Zcentroid, volume, JdensityX, JdensityY, JdensityZ 

Jdensity is optional (the file may be 4-width)
"""
function loadmesh(fn::String; mu_r=1.0)

    data = readdlm(fn, ',', Float64, header=true)[1]

    if size(data)[2] == 4 
        # No Jdensity
        mesh = Mesh(data[:,1:3], data[:,4], mu_r)

    else
        mesh = Mesh(data[:,1:3], data[:,4], data[:,5:7], mu_r)
    end 

    return mesh 
end


"""
    savemesh(fn::String, mesh::Mesh)

Save a finite element mesh to a file. 
"""
function savemesh(fn::String, mesh::Mesh)

    header = ["X [m]" "Y [m]" "Z [m]" "Volume m^3" "JdensityX A/m^2" "JdensityY [A/m^2" "JdensityZ [A/m^3]"]
    N = size(mesh.nodes)[1] 
    data = zeros(N,7) 
    data[:,1:3] = mesh.nodes 
    data[:,4] = mesh.volumes 
    data[:,5:7] = mesh.Jdensity 

    data = vcat(header, data) 

    writedlm(fn,data,",")

end


"""
    loadrings(fn::String)

Load Ring objects from file.

# Returns 
`Vector{CircularRing}, Vector{RectangularRing}`
"""
function loadrings(fn::String)

    data = readdlm(fn, ',', header=true)[1]

    # Read data line by line 
    # This is for flexibility: one file type can load rings of either 
    # rectangular or circular cross-section

    # use parse(Float64, string)

    circ = Vector{CircularRing}(undef,0)
    rect = Vector{RectangularRing}(undef,0)

    for i in range(1, size(data)[1])

        if data[i,4] != "" 
            # This is a circular cross-section 
            push!(circ, CircularRing(string(data[i,1]), data[i,2], data[i,3], data[i,4], data[i,7]))

        else 
            # Rectangular cross-section
            println(i)
            push!(rect, RectangularRing(string(data[i,1]), data[i,2], data[i,3], data[i,5], data[i,6], data[i,7]))

        end
    end

    return circ, rect

end


"""
    saverings(fn::String, rings::Vector{<:Source})

Save Ring objects to a file
"""
function saverings(fn::String, rings::Vector{<:Source})

    header = ["Name" "Z height [m]" "Major Radius [m]" "Minor Radius [m]" "Width [m]" "Height [m]" "Current [A]"]
    data = Matrix{Any}(undef, length(rings), 7)

    for i in eachindex(rings)
        data[i,1] = rings[i].name 
        data[i,3] = rings[i].H 
        data[i,4] = rings[i].R
        data[i,8] = rings[i].I

        if isa(rings[i], CircularRing) 
            data[i,5] = rings[i].r
        else 
            data[i,6] = rings[i].w 
            data[i,7] = rings[i].h
        end 

    end

    data = vcat(header, data) 

    writedlm(fn, data, ',')

end


"""
    loadwires(fn::String)

Load Wire objects from a file
"""
function loadwires(fn::String)

    data = readdlm(fn, ',', header=true)[1]
    N = size(data)[1]
    wires = Vector{Wire}(undef, N) 

    for i in range(1, N) 
        wires[i] = Wire(data[i,1:3], data[i,4:6], data[i,7], data[i,8])
    end

    return wires
end


"""
    savewires(fn::String, wires::Vector{Wire})

Save Wire objects to a file 
"""
function savewires(fn::String, wires::Vector{Wire})

    N = length(wires)
    data = zeros(N, 8)
    header = ["a0_x [m]" "a0_y [m]" "a0_z [m]" "a1_x [m]" "a1_y [m]" "a1_z [m]" "Radius [m]" "Current [A]"]

    for i in eachindex(wires) 
        data[i,1:3] = wires.a0 
        data[i,4:6] = wires.a1 
        data[i,7] = wires.R 
        data[i,8] = wires.I
    end

    writedlm(fn, vcat(header, data), ',')
end
