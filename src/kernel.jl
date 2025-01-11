""" C Kernel for Wired.jl
"""

struct CWire32
    a0::NTuple{3, Cfloat}
    a1::NTuple{3, Cfloat}
    I::Cfloat 
    R::Cfloat
end

struct CWire64
    a0::NTuple{3, Cdouble}
    a1::NTuple{3, Cdouble}
    I::Cdouble 
    R::Cdouble
end

struct CRing32 
	H::Cfloat 
	R::Cfloat 
	r::Cfloat 
	I::Cfloat 
end

struct CRing64 
	H::Cdouble 
	R::Cdouble 
	r::Cdouble 
	I::Cdouble 
end

wires_sp = string(@__DIR__)*"/kernel/"*"wires_sp.so"
wires_dp = string(@__DIR__)*"/kernel/"*"wires_dp.so"
rings_sp = string(@__DIR__)*"/kernel/"*"rings_sp.so"
rings_dp = string(@__DIR__)*"/kernel/"*"rings_dp.so"

""" 
	installkernel()

Install the C kernel by compiling using gcc/make commands
"""
function installkernel()
	current_directory = @__DIR__
	cd(current_directory*"/kernel")
	run(`make`);
	cd(current_directory)
end


"""
	convertCWires(wires::Vector{Wire{Float32}})

Convert Wire objects to CWire objects.
"""
function convertCWires(wires::AbstractArray{Wire{Float32}})
	N = length(wires)
	cwires = Vector{CWire32}(undef, N)
	for i in 1:N 
		cwires[i] = CWire32(
							(wires[i].a0[1], wires[i].a0[2], wires[i].a0[3]), 
							(wires[i].a1[1], wires[i].a1[2], wires[i].a1[3]), 
							wires[i].I, 
							wires[i].R
							)
	end
	
	return cwires 		
end

"""
	convertCWires(wires::Vector{Wire{Float64}})

Convert Wire objects to CWire objects.
"""
function convertCWires(wires::AbstractArray{Wire{Float64}})
	N = length(wires)
	cwires = Vector{CWire64}(undef, N)
	for i in 1:N 
		cwires[i] = CWire64(
							(wires[i].a0[1], wires[i].a0[2], wires[i].a0[3]), 
							(wires[i].a1[1], wires[i].a1[2], wires[i].a1[3]), 
							wires[i].I, 
							wires[i].R
							)
	end
	
	return cwires 		
end

"""
	convertCRingd(wires::Vector{CircularRing{Float32}})

Convert Ring objects to CRing objects.
"""
function convertCRings(rings::AbstractArray{CircularRing{Float32}})
	N = length(rings)
	crings = Vector{CRing32}(undef, N)
	for i in 1:N 
		crings[i] = CRing32(rings[i].H, rings[i].R, rings[i].r, rings[i].I)
	end
	
	return crings
end

"""
	convertCRings(wires::Vector{CircularRing{Float64}})

Convert Ring objects to CRing objects.
"""
function convertCRings(rings::AbstractArray{CircularRing{Float64}})
	N = length(rings)
	crings = Vector{CRing64}(undef, N)
	for i in 1:N 
		crings[i] = CRing64(rings[i].H, rings[i].R, rings[i].r, rings[i].I)
	end
	
	return crings
end

"""
	bs_cwire(nodes::AbstractArray{Float32}, wires::Vector{Wire{Float32}};
					mu_r=1.0)
"""
function bs_cwires(nodes::AbstractArray{Float32}, wires::AbstractArray{Wire{Float32}};
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
	wire_ptr = Base.unsafe_convert(Ptr{CWire32}, cwires)

	# Convert bool to C int
	if check_inside 
		check = 1.0f0 
	else
		check = 0.0f0
	end

	@ccall wires_sp.bfield_wires(Bx_ptr::Ptr{Float32}, 
								   By_ptr::Ptr{Float32}, 
								   Bz_ptr::Ptr{Float32}, 
								   x_ptr::Ptr{Float32},
								   y_ptr::Ptr{Float32},
								   z_ptr::Ptr{Float32}, 
								   wire_ptr::Ptr{CWire32},
								   Nn::Int32, 
								   Nw::Int32, 
								   mu_r::Float32, 
								   check::Int32)::Cvoid
	
	# Zero out singularity points
	map!(x -> isnan(x) ? 0.0 : x, Bx, Bx)
	map!(x -> isnan(x) ? 0.0 : x, By, By)
	map!(x -> isnan(x) ? 0.0 : x, Bz, Bz)
	return hcat(Bx, By, Bz)
end 

"""
	bs_cwire(nodes::AbstractArray{Float64}, wires::Vector{Wire{Float64}};
					mu_r=1.0)
"""
function bs_cwires(nodes::AbstractArray{Float64}, wires::AbstractArray{Wire{Float64}};
					mu_r=1.0)

	Nn = convert(Int64, size(nodes)[1])
	Nw = convert(Int64, length(wires))
	Bx = zeros(Float64, Nn)
	By = zeros(Float64, Nn)
	Bz = zeros(Float64, Nn)
	mu_r = convert(Float64, mu_r)
	cwires = convertCWires(wires)

	Bx_ptr = Base.unsafe_convert(Ptr{Float64}, Bx)
	By_ptr = Base.unsafe_convert(Ptr{Float64}, By)
	Bz_ptr = Base.unsafe_convert(Ptr{Float64}, Bz)
	x_ptr = Base.unsafe_convert(Ptr{Float64}, @view nodes[:,1])
	y_ptr = Base.unsafe_convert(Ptr{Float64}, @view nodes[:,2])
	z_ptr = Base.unsafe_convert(Ptr{Float64}, @view nodes[:,3])
	wire_ptr = Base.unsafe_convert(Ptr{CWire64}, cwires)

	# Convert bool to C int
	if check_inside 
		check = 1.0f0 
	else
		check = 0.0f0
	end

	@ccall wires_dp.bfield_wires(Bx_ptr::Ptr{Float64}, 
								   By_ptr::Ptr{Float64}, 
								   Bz_ptr::Ptr{Float64}, 
								   x_ptr::Ptr{Float64},
								   y_ptr::Ptr{Float64},
								   z_ptr::Ptr{Float64}, 
								   wire_ptr::Ptr{CWire64},
								   Nn::Int64, 
								   Nw::Int64, 
								   mu_r::Float64, 
								   check::Int64)::Cvoid
	
	# Zero out singularity points
	map!(x -> isnan(x) ? 0.0 : x, Bx, Bx)
	map!(x -> isnan(x) ? 0.0 : x, By, By)
	map!(x -> isnan(x) ? 0.0 : x, Bz, Bz)
	return hcat(Bx, By, Bz)
end 


"""
	bs_rings!(nodes::AbstractArray{Float32}, wires::Vector{Wire{Float32}};
					mu_r=1.0)
"""
function bs_crings(nodes::AbstractArray{Float32}, rings::AbstractArray{CircularRing{Float32}};
					mu_r=1.0)

	Nn = convert(Int32, size(nodes)[1])
	Nr = convert(Int32, length(rings))
	Bx = zeros(Float32, Nn)
	By = zeros(Float32, Nn)
	Bz = zeros(Float32, Nn)
	mu_r = convert(Float32, mu_r)
	crings = convertCRings(rings)

	Bx_ptr = Base.unsafe_convert(Ptr{Float32}, Bx)
	By_ptr = Base.unsafe_convert(Ptr{Float32}, By)
	Bz_ptr = Base.unsafe_convert(Ptr{Float32}, Bz)
	x_ptr = Base.unsafe_convert(Ptr{Float32}, @view nodes[:,1])
	y_ptr = Base.unsafe_convert(Ptr{Float32}, @view nodes[:,2])
	z_ptr = Base.unsafe_convert(Ptr{Float32}, @view nodes[:,3])
	ring_ptr = Base.unsafe_convert(Ptr{CRing32}, crings)

	# Convert bool to C int
	if check_inside 
		check = 1.0f0 
	else
		check = 0.0f0
	end

	@ccall rings_sp.bfield_rings(Bx_ptr::Ptr{Float32}, 
								   By_ptr::Ptr{Float32}, 
								   Bz_ptr::Ptr{Float32}, 
								   x_ptr::Ptr{Float32},
								   y_ptr::Ptr{Float32},
								   z_ptr::Ptr{Float32}, 
								   ring_ptr::Ptr{CRing32},
								   Nn::Int32, 
								   Nr::Int32, 
								   mu_r::Float32, 
								   check::Int32)::Cvoid
	
	# Zero out singularity points
	map!(x -> isnan(x) ? 0.0 : x, Bx, Bx)
	map!(x -> isnan(x) ? 0.0 : x, By, By)
	map!(x -> isnan(x) ? 0.0 : x, Bz, Bz)
	# return hcat(Bx, By, Bz)
end 


"""
	bs_rings!(nodes::AbstractArray{Float64}, rings::Vector{Wire{Float64}};
					mu_r=1.0)
"""
function bs_crings(nodes::AbstractArray{Float64}, rings::AbstractArray{CircularRing{Float64}};
					mu_r=1.0)

	Nn = convert(Int64, size(nodes)[1])
	Nr = convert(Int64, length(rings))
	Bx = zeros(Float64, Nn)
	By = zeros(Float64, Nn)
	Bz = zeros(Float64, Nn)
	mu_r = convert(Float64, mu_r)
	crings = convertCRings(rings)

	Bx_ptr = Base.unsafe_convert(Ptr{Float64}, Bx)
	By_ptr = Base.unsafe_convert(Ptr{Float64}, By)
	Bz_ptr = Base.unsafe_convert(Ptr{Float64}, Bz)
	x_ptr = Base.unsafe_convert(Ptr{Float64}, @view nodes[:,1])
	y_ptr = Base.unsafe_convert(Ptr{Float64}, @view nodes[:,2])
	z_ptr = Base.unsafe_convert(Ptr{Float64}, @view nodes[:,3])
	ring_ptr = Base.unsafe_convert(Ptr{CRing64}, crings)

	# Convert bool to C int
	if check_inside 
		check = 1.0f0 
	else
		check = 0.0f0
	end

	@ccall rings_sp.bfield_rings(Bx_ptr::Ptr{Float64}, 
								   By_ptr::Ptr{Float64}, 
								   Bz_ptr::Ptr{Float64}, 
								   x_ptr::Ptr{Float64},
								   y_ptr::Ptr{Float64},
								   z_ptr::Ptr{Float64}, 
								   ring_ptr::Ptr{CRing64},
								   Nn::Int64, 
								   Nr::Int64, 
								   mu_r::Float64, 
								   check::Int64)::Cvoid
	
	# Zero out singularity points
	map!(x -> isnan(x) ? 0.0 : x, Bx, Bx)
	map!(x -> isnan(x) ? 0.0 : x, By, By)
	map!(x -> isnan(x) ? 0.0 : x, Bz, Bz)
	# return hcat(Bx, By, Bz)
end 