""" C Kernel for Wired.jl
"""

# using Libdl 

# kernelpath = string(@__DIR__)*"/kernel/"
# if "ckernel" in readdir(kernelpath)
# 	ckernel = Libdl.dlopen(kernelpath*"ckernel.so")	
# end
kernelpath = string(@__DIR__)*"/kernel/"*"ckernel.so"
""" 
	installkernel()

Install the kernel by compiling using gcc/make commands
"""
function installkernel!(ckernel=Wired.ckernel)

	current_directory = @__DIR__
	cd(current_directory*"/kernel")
	run(`make`);
	cd(current_directory)
	ckernel = Libdl.dlopen(kernelpath*"ckernel.so")
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
	biotsavart_ckernel(nodes::AbstractArray{Float32}, wires::Vector{Wire{Float32}};
					mu_r=1.0)
"""
function biotsavart_ckernel(nodes::AbstractArray{Float32}, wires::Vector{Wire{Float32}};
					mu_r=1.0)

	Nn = convert(Int32, size(nodes)[1])
	Nw = convert(Int32, length(wires))
	B = zeros(Float32, Nn, 3)
	mu_r = convert(Float32, mu_r)

	Bx_ptr = Base.unsafe_convert(Ptr{Float32}, @view B[:,1])
	By_ptr = Base.unsafe_convert(Ptr{Float32}, @view B[:,2])
	Bz_ptr = Base.unsafe_convert(Ptr{Float32}, @view B[:,3])
	x_ptr = Base.unsafe_convert(Ptr{Float32}, @view nodes[:,1])
	y_ptr = Base.unsafe_convert(Ptr{Float32}, @view nodes[:,2])
	z_ptr = Base.unsafe_convert(Ptr{Float32}, @view nodes[:,3])

	@ccall kernelpath.bfield_wires(Bx_ptr::Ptr{Float32}, 
								   By_ptr::Ptr{Float32}, 
								   Bz_ptr::Ptr{Float32}, 
								   x_ptr::Ptr{Float32},
								   y_ptr::Ptr{Float32},
								   z_ptr::Ptr{Float32}, 
								   Ref(wires)::Ref{Vector{Wire{Float32}}},
								   Nn::Int32, 
								   Nw::Int32, 
								   mu_r::Float32)::Cvoid

	return B
end 