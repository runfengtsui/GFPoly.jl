# Euclidean Algorithm
function Base.gcd(f::GFPoly{p1}, g::GFPoly{p2}) where {p1, p2}
    p1 == p2 || throw(ArgumentError("The two polynomials belong to different Polynomial Ring!"))

    r0, r1 = f, g
    iter = 1
    itermax = deg(r1) + 1
    while !iszero(r1) && iter <= itermax
        r2 = rem(r0, r1)
        r0, r1 = r1, r2
        iter += 1
    end
    return r0
end

function isirreducible(f::GFPoly{p}) where p
    n = deg(f)  # the degree of f
    for i = 1:floor(Int, n / 2)
        g_coeff = zeros(GF{p}, p^i+1)
        g_coeff[2], g_coeff[end] = GF{p}(-1), GF{p}(1)
        g = GFPoly{p}(g_coeff)
        r = gcd(g, f)
        deg(r) == 0 || return false
    end
    return true
end

# determine whether a polynomial is primitive
function isprimitive(f::GFPoly{p}) where p
    # first, the polynomial must be irreducible
    isirreducible(f) || return false
    # x^{(q^n-1)/t} - 1 can not be divided by f(x)
    # where t is all factors of p^n - 1
    n = deg(f)  # the degree of f
    factors = factor(Set, p^n - 1)
    for t in factors
        g_coeff = zeros(GF{p}, div(p^n - 1, t)+1)
        g_coeff[1], g_coeff[end] = GF{p}(-1), GF{p}(1)
        g = GFPoly{p}(g_coeff)
        iszero(rem(g, f)) && return false
    end
    return true
end

