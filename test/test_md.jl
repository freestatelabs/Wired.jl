

module TestMD
export Test, f
precision = Float32

struct Test{T<:Real}

	a::T
	b::AbstractArray{T}

	function Test{T}(a::Real, b::Vector{<:Real}) where T

		new{T}(convert(T, a), convert.(T,b))
	end

	function Test(a::Real, b::Vector{<:Real})

		new{precision}(convert(precision, a), convert.(precision, b))
	end

end


function f(t::Vector{Test{T}}) where T<:Real
	println(t[1].a)
	println(t[1].b)
end

t = Test(1, [2,3])
f([t])
t2 = Test{Int8}(1, [2,3])
f([t2])
end