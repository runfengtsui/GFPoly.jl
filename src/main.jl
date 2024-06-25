using Parsers

include("types.jl")
include("arithmetic.jl")
include("algorithm.jl")

# Brute force search to find all n-th degree irreducible polynomials in GF(p)
function irrepoly(p::T, n::T) where T <: Integer
    poly = GFPoly{p}[]
    Fp = GF{p}
    for k = p^n:2*p^n-1
        padic = string(k, base=p)
        f_coeff = reverse(parse.(Int, split(padic, "")))
        f_coeff = map(Fp, f_coeff)
        f = GFPoly{p}(f_coeff)
        isirreducible(f) && push!(poly, f)
    end
    return poly
end

# Brute force search to find all n-th degree primitive polynomials in GF(p)
function primpoly(p::T, n::T) where T <: Integer
    poly = GFPoly{p}[]
    for k = p^n:2*p^n-1
        padic = string(k, base=p)
        f_coeff = reverse(parse.(Int, split(padic, "")))
        f = GFPoly{p}(map(GF{p}, f_coeff))
        isprimitive(f) && push!(poly, f)
    end
    return poly
end

