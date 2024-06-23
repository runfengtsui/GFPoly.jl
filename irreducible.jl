using Polynomials
using Parsers

# call C libgmp's __gmpz_probab_prime_p function to test prime
function is_probable_prime(x::Integer, reps::Integer=25)
    return ccall((:__gmpz_probab_prime_p, :libgmp), Cint,
                 (Ref{BigInt}, Cint), x, reps) != 0
end

# 判断两个多项式在模 p 情况下是否整除
function isdivided(coef::Vector{T}, coeg::Vector{T}, p::T) where T <: Integer
    f = Polynomial(coef)
    g = Polynomial(coeg)
    r = rem(f, g)
    all(rem.(r.coeffs, p) .== 0) && return true
    return false
end

# 检测是否有一次因式
function isdividedbyone(coef::Vector{T}, p::T) where T <: Integer
    for a = 0:p-1
        # 一次多项式的系数
        coeg = [a, 1]
        isdivided(coef, coeg, p) && return true
    end
    return false
end

function get_n_irrepoly_Fp(p::T, n::T) where T <: Integer
    # 有限域 GF(p^n) 中 p 必须是素数
    !is_probable_prime(p) && throw(DomainError(p, "Characteristic is not prime in GF(p)"))

    # 保存所有 n 次不可约多项式
    irrepoly = Vector{Int64}[]
    # construct x^(p^n) - x
    coef = zeros(T, p^n+1)
    coef[2] = -1 # -x
    coef[end] = 1 # 首一多项式, 所有不可约多项式的乘积

    for k = p^n:2*p^n-1
        # n 次首一多项式的系数的 p 进制表示
        padic = string(k, base=p)
        # 将 p 进制转化为系数向量
        coeg = reverse(parse.(Int64, split(padic, "")))
        if n == 1
            # 一次多项式都是不可约多项式
            push!(irrepoly, coeg)
            continue
        end
        # 待检测 n 次多项式是否存在一次因式
        isdividedbyone(coeg, p) && continue
        # 待检测多项式是否可以整除 x^(p^n) - x
        isdivided(coef, coeg, p) && push!(irrepoly, coeg)
    end
    return irrepoly
end

