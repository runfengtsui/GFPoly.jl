using Polynomials
using Parsers
using Primes

# 素因子分解
function prime_factors(n::T) where T <: Integer
    factors = T[]
    d = BigInt(2)
    while n > 1
        # 如果 d 是因子，则从 n 中都除掉
        while n > 1 && n % d == 0 
            n = div(n, d)
            push!(factors, d)
        end

        if n == 1 # 已经没有素因子
            break
        else
            # 考虑下一个可能的因子，所有因子都 > d了
            d += (d == 2 ? 1 : 2)
            if n < d^2
                # 如果n还是合数，n=ab, 则a>=d, b>=d, n >= d^2
                # 所以n < d^2时n为素数
                push!(factors, n)
                break
            end
            # 这里n还有可能是合数，继续用d尝试能否除尽
        end # if - else
    end # while n > 1
    return factors
end

# 判断两个多项式在模 p 情况下是否整除
function isdivided(f::Polynomial{T}, g::Polynomial{T}, p::T) where T <: Integer
    r = rem(f, g)
    all(rem.(r.coeffs, p) .== 0) && return true
    return false
end

# for all ni in {n/p1, n/p2, ..., n/pt}, gcd(x^{q^ni-1}-1, f(x)) = 1
function isgcdeqone(f::Polynomial{T}, p::T, n::T) where T <: Integer
    # M = div.(n, Set(prime_factors(n))) # 使用自己的素因子分解函数
    M = div.(n, factor(Set, n))
    for ni in M
        coeh = zeros(T, p^(ni-1))
        coeh[1], coeh[end] = -1, 1
        h = Polynomial(coeh)
        gcdpoly = gcd(f, h)
        # 系数向量的长度为 1, 说明 gcd(g, h) = 1
        length(gcdpoly.coeffs) != 1 && return false
    end
    return true
end

function isirreducible(f::Polynomial{T}, p::T) where T <: Integer
    # 获取多项式的次数
    n = degree(f)
    # 所有一次多项式都是不可约多项式
    n == 1 && return true
    # 构造首一多项式, 所有不可约多项式的乘积 x^(p^n) - x 
    coeg = zeros(T, p^n+1)
    coeg[2], coeg[end] = -1, 1
    g = Polynomial(coeg)
    # condition 1: f(x) | x^(p^n) - x
    # condition 2: for all c in GF(P), f(c) != 0
    # condition 3: for all ni in {n/p1, n/p2, ..., n/pt}, gcd(x^{q^ni-1}-1, f(x)) = 1
    isdivided(g, f, p) && all(rem.(f.(0:p-1), p) .!= 0) && isgcdeqone(f, p, n) && return true
    return false
end

function get_n_irrepoly_Fp(p::T, n::T) where T <: Integer
    # 有限域 GF(p^n) 中 p 必须是素数
    !isprime(p) && throw(DomainError(p, "Characteristic is not prime in GF(p)"))
    # 保存所有 n 次不可约多项式
    irrepoly = Polynomial{T}[]

    for k = p^n:2*p^n-1
        # n 次首一多项式的系数的 p 进制表示
        padic = string(k, base=p)
        # 将 p 进制转化为系数向量
        coeg = reverse(parse.(Int64, split(padic, "")))
        g = Polynomial(coeg)
        println(g)
        isirreducible(g, p) && push!(irrepoly, g)
    end
    return irrepoly
end

# 验证是否为本原多项式比不可约多项式多的一个条件
# f(x) 不能整除 x^{(q^n-1)/t} - 1
function isprimitive(f::Polynomial{T}, p::T, n::T) where T <: Integer
    factors = factor(Set, p^n - 1)
    for t in factors
        deg = div(p^n - 1, t)
        coeg = zeros(T, deg)
        coeg[1], coeg[end] = -1, 1
        g = Polynomial(coeg)
        isdivided(g, f, p) && return false
    end
    return true
end

function get_n_pripoly_Fp(p::T, n::T) where T <: Integer
    # 有限域 GF(p^n) 中 p 必须是素数
    !isprime(p) && throw(DomainError(p, "Characteristic is not prime in GF(p)"))

    # 保存所有 n 次不可约多项式
    pripoly = Polynomial{T}[]
    coef = zeros(T, p^n+1)
    coef[2], coef[end] = -1, 1
    f = Polynomial(coef)

    for k = p^n:2*p^n-1
        padic = string(k, base=p)
        coeg = reverse(parse.(Int64, split(padic, "")))
        g = Polynomial(coeg)
        if n == 1
            push!(pripoly, g)
            continue
        end
        # g(x) | x^(p^n) - x
        !isdivided(f, g, p) && continue
        # for all c in GF(P), g(c) != 0
        any(rem.(g.(0:p-1), p) .== 0) && continue
        # for all ni in {n/p1, n/p2, ..., n/pt}, gcd(x^{q^ni-1}-1, f(x)) = 1
        !isgcdeqone(g, p, n) && continue
        isprimitive(g, p, n) && push!(pripoly, g)
    end
    return pripoly
end

