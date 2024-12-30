""" Utility functions for Wired.jl
    (c) 2024 ryan@freestatelabs
"""


"""
    function ellipK(k2; itmax=100, errmax=1e-8)

Solve the complete elliptic integral of the first kind
"""
function ellipK(k2::Real; itmax=100, errmax=1e-8)

    if abs(k2 - 1.0) <= errmax
        return Inf 
    end

    it = 1 
    err = 2*errmax

    a0 = 1
    g0 = sqrt(1-k2)

    while it < itmax && err > errmax

        a1 = 0.5*(a0+g0)
        g1 = sqrt(a0*g0)

        a0 = a1 
        g0 = g1

        err = abs(a0-g0)
        it += 1

    end

    return pi/(2*a0)

end

"""
    function ellipE(k2; itmax=100, errmax=1e-8)

Solve the complete elliptic integral of the second kind

Reference:
[https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html](https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html)
"""
function ellipE(k2::Real; itmax=100, errmax=1e-8)

    E = 0.0

    if abs(k2 - 1.0) <= errmax
        E = 1.0 
    elseif abs(k2 - 0.0) <= errmax 
        E = pi/2 
    else 

        err = 2*errmax
        n = 0
        an = 1 
        gn = sqrt(1-k2)
        cn = abs(an^2 - gn^2)
        esumn = cn*2.0^(n-1)

        while n < itmax && err > errmax
            n +=1 

            an1 = (an + gn)/2.0
            gn1 = sqrt(an*gn)
            cn1 = abs(an1^2 - gn1^2)
            esumn1 = esumn + cn1*(2^(n-1))

            err = abs(esumn1 - esumn)

            an = an1 
            gn = gn1 
            esumn = esumn1 

        end

        E = (1-esumn)*pi/(2*an)
    end

    return E

end


"""
    ellipKE(k2; itmax=100, errmax=1e-12)

Solve the complete elliptic integrals of the first and second kinds.
"""
function ellipKE(k2::Real; itmax=100, errmax=1e-12)

    k2 = Float64(k2)

    if abs(k2 - 1.0) < errmax 
        K = Inf 
        E = 1.0 

    elseif abs(k2) < errmax 
        K = pi/2 
        E = pi/2 

    else
        n = 0
        err = 2*errmax
        
        a0 = 1 
        g0 = sqrt(1-k2)
        c0 = abs(a0^2 - g0^2)
        esum0 = (2.0^(n-1))*c0

        while n < itmax && err > errmax

            a1 = 0.5*(a0+g0)
            g1 = sqrt(a0*g0)
            c1 = sqrt(abs(a1^2 - g1^2))
            esum1 = esum0 + (2.0^(n-1))*c1^2
    
            err = max(abs(a1-g1),abs(esum1-esum0))
            a0 = a1 
            g0 = g1
            esum0 = esum1
            n += 1
    
        end

        K = pi/(2*a0)
        E = (pi/(2*a0))*(1-esum0)

    end

    return K, E

end


"""
    crossrows!(C::AbstractArray, A::AbstractArray, b::Vector)

Cross each row of A by a 3-length vector b and place in C, i.e. 

    C[i,:] = cross(A[i,:], b)

"""
@views function crossrows!(C::AbstractArray, A::AbstractArray, b::AbstractVector)

    C[:,1] .= A[:,2] .* b[3] .- A[:,3] .* b[2]
    C[:,2] .= A[:,3] .* b[1] .- A[:,1] .* b[3]
    C[:,3] .= A[:,1] .* b[2] .- A[:,2] .* b[1]

end


"""
    normrows!(b::Vector, A::AbstractArray)

Calculate the vector norm of each row in Nx3 matrix and place in `b`, i.e. 

    b[i] = norm(A[i,:])

"""
@views function normrows!(b::AbstractVector, A::AbstractArray)

    b .= (sqrt.(A[:,1].^2 .+ A[:,2].^2 .+ A[:,3].^2))

end


"""
    multrows!(A::AbstractArray, b::Vector)

In-place element-wise multiplication of each row of an Nx3 matrix by the 
corresponding element in vector `b`, i.e.

    A[i,:] .*= b[i]

"""
@views function multrows!(A::AbstractArray, b::AbstractVector)

    A[:,1] .*= b 
    A[:,2] .*= b 
    A[:,3] .*= b

end


""" 
    dotrows!(c::Vector, A::AbstractArray, b::Vector)

In-place dot-product of the rows of A by vector b; place in C, i.e. 

    c[i] .*= dot(A[i,:], b)       # or C[i,:] .*= sum(A[i,:] .* b)
"""
@views function dotrows!(c::AbstractVector, A::AbstractArray, b::AbstractVector)

    c .= (A[:,1] .* b[1] .+ A[:,2] .* b[2] .+ A[:,3] .* b[3])
end


# from: https://stackoverflow.com/questions/51466537/how-to-plot-a-vector-field-in-julia
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))


"""
    threadindices(it::Integer, Nt::Integer, N::Integer) 

Split a problem up for multi-threaded operation

For thread number `it`, total number of threads `Nt`, and total number of 
tasks `N`, determine the start/stop index for that particular thread number
"""
function threadindices(it::Integer, Nt::Integer, N::Integer) 

    Nperthread = div(N, Nt)
    remainder = rem(N, Nt) 

    if it == 1
        i1 = 1 
        i2 = i1 + Nperthread - 1 
    elseif it == Nt 
        i2 = N 
        i1 = i2 - Nperthread - remainder + 1 
    else
        i1 = (it-1)*Nperthread + 1
        i2 = i1 + Nperthread - 1 
    end 

    return i1:i2
end


"""
    matrixtotable(A::AbstractArray, header::AbstractArray, digits=1)

    Convert a `Matrix` into a Markdown table

    # Returns
    `String` that can be printed and copied into a markdown document
"""
function matrixtotable(A::AbstractArray, header::AbstractArray, digits=1)

    output = ""

    for h in header
        output *= "| " * h * " "
    end 
    output *= "|\n"

    for i in 1:size(header)[1]
        output *= "| --- "
    end 
    output *= " |\n"

    for i in 1:size(A)[1] 
        for j in 1:size(A)[2]

            output *= "| " * string(round(A[i,j], digits=digits)) * " "

        end
        output *= " |\n"
    end

    return output

end