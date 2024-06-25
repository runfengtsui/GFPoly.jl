import Base: +, -, *, /, ==, <

# Arithmetic operations over GF(p)
# add two elements in the same Galois Field GF(p)
function +(a::GF{p1}, b::GF{p2})::GF{p1} where {p1, p2}
    p1 == p2 || throw(ArgumentError("Galois Field have different characteristic!"))
    m = p1 - b.x
    s = a.x < m ? a.x + b.x : a.x - m
    return GF{p1}(s)
end

# get the inverse element of addition
-(a::GF{p}) where p = a.x == 0 ? deepcopy(a) : GF{p}(p - a.x)
# minus two elements in the same Galois Field GF(p)
function -(a::GF{p1}, b::GF{p2})::GF{p1} where {p1, p2}
    p1 == p2 || throw(ArgumentError("Galois Field have different characteristic!"))
    m = a.x < b.x ? a.x + (p1 - b.x) : a.x - b.x
    return GF{p1}(m)
end

# calculate the multiplication of two elements in the same Galois Field GF(p)
function *(a::GF{p1}, b::GF{p2})::GF{p1} where {p1, p2}
    p1 == p2 || throw(ArgumentError("Galois Field have different characteristic!"))
    return GF{p1}(widemul(a.x, b.x))
end

# get the inverse element of multiplication over GF*(p)
function Base.inv(a::GF{p}) where p
    a.x == 0 && throw(DivideError())
    g, s, _ = gcdx(a.x, p)
    g != 1 && error("Characteristic not prime in GF_", p)
    return GF{p}(s)
end

# calculate the division of two elements in the same Galois Field GF(p)
function /(a::GF{p1}, b::GF{p2})::GF{p1} where {p1, p2}
    p1 == p2 || throw(ArgumentError("Galois Field have different characteristic!"))
    return a * inv(b)
end

# isequal
==(a::GF{p1}, b::GF{p2}) where {p1, p2} = (p1 == p2) && (a.x == b.x)

# is less than
function <(a::GF{p1}, b::GF{p2}) where {p1, p2}
    p1 == p2 || throw(ArgumentError("Galois Field have different characteristic!"))
    return a.x < b.x ? true : false
end

###########################################################################################
# Arithmetic operations over GF(p)[x]

# calculate the addition of two polynomials over GF(p)[x]
function +(f::GFPoly{p1}, g::GFPoly{p2}) where {p1, p2}
    p1 == p2 || throw(ArgumentError("The two polynomials belong to different Polynomial Ring!"))
    f_coeff, g_coeff = f.coeff, g.coeff
    n = max(length(f_coeff), length(g_coeff))
    c = [get(f_coeff, i, 0) + get(g_coeff, i, 0) for i in 1:n]
    while n > 1 && c[n] == 0
        pop!(c)
        n -= 1
    end
    return GFPoly{p1}(c)
end

# get the inverse element of addition over GF(p)[x]
-(f::GFPoly{p}) where p = GFPoly{p}(-f.coeff)
# calculate the minus of two polynomials over GF(p)[x]
function -(f::GFPoly{p1}, g::GFPoly{p2}) where {p1, p2}
    p1 == p2 || throw(ArgumentError("The two polynomials belong to different Polynomial Ring!"))
    return f + (-g)
end

# number multiplication
*(k::GF{p}, f::GFPoly{p}) where p = GFPoly{p}(k * f.coeff)
*(f::GFPoly{p}, k::GF{p}) where p = GFPoly{p}(k * f.coeff)
*(k::Integer, f::GFPoly{p}) where p = GF{p}(k) * f
*(f::GFPoly{p}, k::Integer) where p = GF{p}(k) * f
# calculate the multiplication of two polynomials over GF(p)[x]
function *(f::GFPoly{p1}, g::GFPoly{p2}) where {p1, p2}
    p1 == p2 || throw(ArgumentError("The two polynomials belong to different Polynomial Ring!"))
    f_deg, g_deg = deg(f), deg(g)
    f_coeff, g_coeff = f.coeff, g.coeff
    c = zeros(GF{p1}, f_deg + g_deg + 1)
    @simd for i = 0:f_deg
        @simd for j = 0:g_deg
            @inbounds c[i+j+1] += f_coeff[i+1] * g_coeff[j+1]
        end
    end
    return GFPoly{p1}(c)
end

# calculate the quotient and divider and remainder of two polynomials when dividing over GF(p)[x]
function Base.divrem(f::GFPoly{p1}, g::GFPoly{p2}) where {p1, p2}
    p1 == p2 || throw(ArgumentError("The two polynomials belong to different Polynomial Ring!"))
    f_deg, g_deg = deg(f), deg(g)
    q_deg = f_deg - g_deg
    q_deg < 0 && return zero(GFPoly{p1}), f

    g_coeff = g.coeff
    t = inv(g_coeff[end])
    q_coeff = zeros(GF{p1}, q_deg+1)
    r_coeff = deepcopy(f.coeff)
    @simd for i = q_deg:-1:0
        @inbounds q_coeff[i+1] = t * r_coeff[i+g_deg+1]
        @simd for j = 0:g_deg
            @inbounds r_coeff[i+j+1] -= q_coeff[i+1] * g_coeff[j+1]
        end
    end
    return GFPoly{p1}(q_coeff), GFPoly{p1}(r_coeff)
end

# get the remainder
Base.rem(f::GFPoly{p1}, g::GFPoly{p2}) where {p1, p2} = divrem(f, g)[2]

