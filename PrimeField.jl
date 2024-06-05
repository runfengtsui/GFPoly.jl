import Base: +, -, *, ==

# call C libgmp's __gmpz_probab_prime_p function to test prime
function is_probable_prime(x::Integer, reps::Integer=25)
    return ccall((:__gmpz_probab_prime_p, :libgmp), Cint,
                 (Ref{BigInt}, Cint), x, reps) != 0
end

# Type: Finite Field Zp where p is a prime
# struct PrimeField{T <: Integer}
#     p::T
#     # constructor for Prime Field Zp
#     function PrimeField(p::T) where T <: Integer
#         !is_probable_prime(p) && throw(DomainError(p, "Characteristic is not prime in PrimeField(p)"))
#         new{T}(p)
#     end
# end
struct PrimeField{T <: Integer}
    p::T
end

# Set the output of PrimeField when printing
Base.show(io::IO, Z::PrimeField) = print(io, "Prime Field Z_", Z.p)

# the characteristic of Zp
char(Z::PrimeField{T}) where T <: Integer = Z.p

# the order of Zp
order(Z::PrimeField{T}) where T <: Integer = Z.p

# interface to construct Prime Field Zp
function GF(p::T) where T <: Integer
    !is_probable_prime(p) && throw(DomainError(p, "Characteristic is not prime in GF(p)"))
    return PrimeField{T}(p)
end

##############################################################################

struct PFElem{T <: Integer}
    e::T
    field::PrimeField{T}
end

# Set the output of GFElem when printing
Base.show(io::IO, x::PFElem) = print(io, x.e)

# Get the value of the element
data(a::PFElem) = a.e

# Get the PFField that the element belongs to
field(a::PFElem) = a.field

function (Z::PrimeField{T})(a::Integer) where T <: Integer
    p = Z.p::T
    d = convert(T, a % p)::T
    # allow a is a negative integer, so d may be negative
    (d < 0) && (d += p)
    return PFElem{T}(d, Z)
end

function (Z::PrimeField{T})(a::PFElem{T}) where T <: Integer
    field(a) != Z && error("The element ", a, " does not belong to the field Z_", Z.p, ".")
    return a
end

##############################################################################

function is_same_field(x::PFElem, y::PFElem)
   x.field != y.field && error("Operations on distinct finite fields not supported")
end

# calculate the inverse element
function -(x::PFElem{T}) where T <: Integer
    if x.e == 0
        return deepcopy(x)
    else
        Z = field(x)
        return PFElem{T}(Z.p - x.e, Z)
    end
end

# add two elements in the same GFField
function +(x::PFElem{T}, y::PFElem{T}) where T <: Integer
    is_same_field(x, y)
    Z = field(x)
    p = char(Z)::T
    m = p - y.e
    s = x.e < m ? x.e + y.e : x.e - m
    return PFElem{T}(s, Z)
end

# minus two elements in the same GFField
function -(x::PFElem{T}, y::PFElem{T}) where T <: Integer
    is_same_field(x, y)
    Z = field(x)
    p = char(Z)::T
    m = x.e < y.e ? x.e + (p - y.e) : x.e - y.e
    return PFElem{T}(m, Z)
end

function *(x::PFElem{T}, y::PFElem{T}) where T <: Integer
    is_same_field(x, y)
    Z = field(x)
    return Z(widemul(x.e, y.e))
end

function *(x::Integer, y::PFElem{T}) where T <: Integer
    Z = field(y)
    return Z(widemul(x, y.e))
end

function *(x::PFElem{T}, y::Integer) where T <: Integer
    Z = field(x)
    return Z(widemul(x.e, y))
end

function ==(x::PFElem{T}, y::PFElem{T}) where T <: Integer
    is_same_field(x, y)
    return x.e == y.e
end

