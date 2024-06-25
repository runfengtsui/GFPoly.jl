using Primes

# Galois Field GF(p) where p is a prime
struct GF{p} <: Integer
    x::UInt64
    # generate an element in Galois Field GF(p)
    function GF{p}(x) where {p}
        !isprime(p) && throw(DomainError(p, "Characteristic is not prime in PrimeField(p)"))
        # allow x is a negative integer, so d may be negative
        d = convert(UInt64, ((x % p) + p) % p)
        new{p}(d)
    end
end
GF{p}(x::GF{p}) where {p} = x

# Set the output of GFElem when printing
Base.show(io::IO, ::Type{GF{p}}) where p = print(io, "GF", p)
Base.show(io::IO, a::GF{p}) where p = print(io, a.x)
Base.string(a::GF{p}) where p = return string(a.x)
Base.promote_rule(::Type{<:Integer}, ::Type{GF{p}}) where p = GF{p}
Base.promote_rule(::Type{GF{p}}, ::Type{<:Integer}) where p = GF{p}
Base.zero(::GF{p}) where p = GF{p}(0)
Base.zero(::Type{GF{p}}) where p = GF{p}(0)
Base.one(::Type{GF{p}}) where p = GF{p}(1)

###########################################################################################

struct GFPoly{p}
    coeff::Vector{GF{p}}
    function GFPoly{p}(coeff) where p
        n = length(coeff)
        while n > 1 && coeff[n] == 0
            pop!(coeff)
            n -= 1
        end
        new{p}(coeff)
    end
end
# Set the output of GFPoly when printing
function Base.show(io::IO, f::GFPoly{p}) where p
    @inbounds for (i, ai) in enumerate(f.coeff)
        if ai > 0
            if i == 1
                print(io, ai)   # print constant
            elseif i == 2
                print(io, " + " , ai == 1 ? "" : string(ai), "x")
            else
                print(io, " + " , ai == 1 ? "" : string(ai), "x^$(i-1)")
            end
        else
            continue
        end
    end
end

# get the degree of the polynomial
deg(f::GFPoly{p}) where p = length(f.coeff) - 1

# zero polynomial
Base.zero(::GFPoly{p}) where p = GFPoly{p}([GF{p}(0)])
Base.zero(::Type{GFPoly{p}}) where p = GFPoly{p}([GF{p}(0)])
# determine whether a polynomial is zero polynomial
Base.iszero(f::GFPoly{p}) where p = deg(f) == 0 && (f.coeff)[1] == GF{p}(0)

# evaluate a polynomial at x
(f::GFPoly{p})(x::GF{p}) where p = Base.Math.evalpoly(x, f.coeff)
(f::GFPoly{p})(x::Integer) where p = f(GF{p}(x))

