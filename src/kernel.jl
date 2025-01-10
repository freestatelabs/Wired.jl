""" C Kernel for Wired.jl
"""

using Libdl

struct CWire
    a0::NTuple{3, Cfloat}
    a1::NTuple{3, Cfloat}
    I::Cfloat 
    R::Cfloat
end

kernelpath = string(@__DIR__)*"/kernel/"*"ckernel.so"

""" 
	installkernel()

Install the kernel by compiling using gcc/make commands
"""
function installkernel()
	current_directory = @__DIR__
	cd(current_directory*"/kernel")
	run(`make`);
	cd(current_directory)
end


"""
	testkernel(a::AbstractArray{Float32})

Test the kernel
"""
function testkernel(a::AbstractArray{Float32})
	b = zeros(Float32, 1)
	println("Value of b before @ccall is: "*string(b))
	@ccall "./kernel/kernel.so".test(a::Ptr{Float32}, b::Ptr{Float32})::Cvoid
	println("Value of b after @ccall is:  "*string(b))
end

"""
	convertCWires(wires::Vector{Wire{Float32}})

Convert Wire objects to CWire objects.
"""
function convertCWires(wires::AbstractArray{Wire{Float32}})
	N = length(wires)
	cwires = Vector{CWire}(undef, N)
	for i in 1:N 
		cwires[i] = CWire(
							(wires[i].a0[1], wires[i].a0[2], wires[i].a0[3]), 
							(wires[i].a1[1], wires[i].a1[2], wires[i].a1[3]), 
							wires[i].I, 
							wires[i].R
							)
	end
	
	return cwires 		
end

"""
	biotsavart_ckernel(nodes::AbstractArray{Float32}, wires::Vector{Wire{Float32}};
					mu_r=1.0)
"""
function biotsavart_ckernel(nodes::AbstractArray{Float32}, wires::AbstractArray{Wire{Float32}};
					mu_r=1.0)

	Nn = convert(Int32, size(nodes)[1])
	Nw = convert(Int32, length(wires))
	Bx = zeros(Float32, Nn)
	By = zeros(Float32, Nn)
	Bz = zeros(Float32, Nn)
	mu_r = convert(Float32, mu_r)
	cwires = convertCWires(wires)

	Bx_ptr = Base.unsafe_convert(Ptr{Float32}, Bx)
	By_ptr = Base.unsafe_convert(Ptr{Float32}, By)
	Bz_ptr = Base.unsafe_convert(Ptr{Float32}, Bz)
	x_ptr = Base.unsafe_convert(Ptr{Float32}, @view nodes[:,1])
	y_ptr = Base.unsafe_convert(Ptr{Float32}, @view nodes[:,2])
	z_ptr = Base.unsafe_convert(Ptr{Float32}, @view nodes[:,3])
	wire_ptr = Base.unsafe_convert(Ptr{CWire}, cwires)

	if check_inside 
		check = 1 
	else
		check = 0 
	end

	@ccall kernelpath.bfield_wires(Bx_ptr::Ptr{Float32}, 
								   By_ptr::Ptr{Float32}, 
								   Bz_ptr::Ptr{Float32}, 
								   x_ptr::Ptr{Float32},
								   y_ptr::Ptr{Float32},
								   z_ptr::Ptr{Float32}, 
								   wire_ptr::Ptr{CWire},
								   Nn::Int32, 
								   Nw::Int32, 
								   mu_r::Float32, 
								   check::Int32)::Cvoid
	map!(x -> isnan(x) ? 0.0 : x, Bx, Bx)
	map!(x -> isnan(x) ? 0.0 : x, By, By)
	map!(x -> isnan(x) ? 0.0 : x, Bz, Bz)
	return hcat(Bx, By, Bz)
end 