using StaticArrays


A = Float32.([1,1,1,1])
B = Float32.([1,2,3,4])
C = Float32.([0,0,0,0])

A_ptr = Base.unsafe_convert(Ptr{Float32}, A)
B_ptr = Base.unsafe_convert(Ptr{Float32}, B)
C_ptr = Base.unsafe_convert(Ptr{Float32}, C)

N = Int32(4)
M = Int32(2)

println(C)
@ccall "./test.so".mmult(A_ptr::Ptr{Float32}, B_ptr::Ptr{Float32}, C_ptr::Ptr{Float32}, N::Int32, M::Int32)::Cvoid
println(C)



abstract type Source end 
precision = Float32

# struct Wire{T<:AbstractFloat} <: Source
#     a0::SVector{3, T}
#     a1::SVector{3, T}
#     I::T
#     R::T

#     # Constructor converts inputs to StaticVectors
#     # Requires specifying a type
#     function Wire{T}(a0::Vector{<:Real}, a1::Vector{<:Real}, I::Real, R::Real) where T<:AbstractFloat
#         new{T}(SVector{3}(convert.(T, a0)), SVector{3}(convert.(T, a1)), convert(T, I), convert(T, R))
#     end

#     # Constructor applies Wired.precision by default (convenience method)
#     function Wire(a0::Vector{<:Real}, a1::Vector{<:Real}, I::Real, R::Real)
#         Wire{precision}(a0, a1, I, R)
#     end
# end


# Define the Wire struct with appropriate layout
struct Wire
    a0::NTuple{3, Cfloat}  # Static tuple of floats
    a1::NTuple{3, Cfloat}  # Static tuple of floats
    I::Cfloat  # Current
    R::Cfloat  # Resistance
end

# Create the Wire struct and an array of Wire structs
wire = Wire((1.0f0, 2.0f0, 3.0f0), (4.0f0, 5.0f0, 6.0f0), 7.0f0, 8.0f0)
wires = [wire, wire, wire]  # Example array of Wire structs

# Allocate a buffer for the result array (8 values per Wire)
a = zeros(Float32, 8 * length(wires))  # Adjust the size based on Nwires
a_ptr = Base.unsafe_convert(Ptr{Float32}, a)

wire_refs = [Ref(wire) for wire in wires]

# Allocate an array of pointers to the references
wire_pointers = Vector{Ptr{Wire}}(undef, length(wire_refs))
for i in 1:length(wire_refs)
    wire_pointers[i] = pointer_from_objref(wire_refs[i])  # Convert Ref{Wire} to Ptr{Wire}
end


# Convert the array of pointers to Ptr{Ptr{Wire}}
wire_ptrs = pointer_from_objref(wire_pointers)

# Call the C function with the appropriate arguments
Nwires = length(wires)
@ccall "./test.so".testwires(wire_ptrs::Ptr{Ptr{Wire}}, Nwires::Int32, a_ptr::Ptr{Float32})::Cvoid

# Print the result
println(a)


println("testing vdot")

a = Float32.([1,2,3,4])
b = Float32.([1,1,1,1])
ab = Float32.([1 1; 2 1; 3 1; 4 1])
c = zeros(Float32, 4,2) 
N = length(wires)
a_ptr = Base.unsafe_convert(Ptr{Float32}, @view ab[:,1])
b_ptr = Base.unsafe_convert(Ptr{Float32}, @view ab[:,2])
c_ptr = Base.unsafe_convert(Ptr{Float32}, @view c[:,1])
display(c)
@ccall "./test.so".vdot(a_ptr::Ptr{Float32}, b_ptr::Ptr{Float32}, c_ptr::Ptr{Float32}, N::Int32)::Cvoid 
display(c)