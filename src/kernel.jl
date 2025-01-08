""" C Kernel for Wired.jl
"""

using Libdl 
println(string(@__DIR__)*"/kernel/kernel.so")
ckernel = Libdl.dlopen(string(@__DIR__)*"/kernel/kernel.so")	
dump(ckernel)

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
	biotsavart_ckernel(nodes::AbstractArray{Float32}, wires::Vector{Wire{Float32}};
					mu_r=1.0)
"""
function biotsavart_ckernel(nodes::AbstractArray{Float32}, sources::AbstractArray{Float32};
					mu_r=1.0)

	dump(ckernel)


	Nn = size(nodes)[1]
	Nw = size(sources)[1]
	B = zeros(Nn, 3)

	B_ptr = Base.unsafe_convert(Ptr{Ptr{Float32}}, B)
	nodes_ptr = Base.unsafe_convert(Ptr{Ptr{Float32}}, nodes)
	sources_ptr = Base.unsafe_convert(Ptr{Ptr{Float32}}, sources)

	@ccall "./kernel/kernel.so".solve2(B_ptr::Ptr{Ptr{Float32}}, nodes_ptr::Ptr{Ptr{Float32}}, sources_ptr::Ptr{Ptr{Float32}}, mu_r::Float32, Nn::Int32, Nw::Int32)::Cvoid
	# ccall((:solve2, "./kernel/kernel.so"), Cvoid, (Ptr{Ptr{Float32}}, Ptr{Ptr{Float32}}, Ptr{Ptr{Float32}}, Float32, Int32, Int32), (B_ptr, nodes_ptr, sources_ptr, mu_r, Nn))
	return B
end 