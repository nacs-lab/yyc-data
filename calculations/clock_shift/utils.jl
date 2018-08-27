#

if VERSION >= v"0.7.0"
    @eval using SpecialFunctions
end

@inline uninit_ary(::Type{T}, args...) where T =
    @static VERSION >= v"0.7.0" ? T(undef, args...) : T(args...)

const _hermit_coeff_cache = Vector{Float64}[[1.0], [1.0]]

# Return coefficient array for the non-zero ones devided by 2^n
# Highest power first.
# If the highest power (n) is even, the return value has (n ÷ 2 + 1) terms
# corresponding to the n-th, (n-2)-th, ... 0-th power.
# If the highest power (n) is odd, the return value has ((n - 1) ÷ 2 + 1) terms
# corresponding to the n-th, (n-2)-th, ... 1st power.
function hermite_coeff(n)
    if n < length(_hermit_coeff_cache)
        return _hermit_coeff_cache[n + 1]
    end
    # H_n+1(x) = 2x H_n(x) - 2n H_n-1(x)
    # H_n+1(x)/2^(n+1) = x H_n(x)/2^(n) - n H_n-1(x)/2^n
    # H'_n+1(x) = x H'_n(x) - n H'_n-1(x)/2
    # H'_n(x) = x H'_n-1(x) - (n - 1) H'_n-2(x)/2
    c_1 = hermite_coeff(n - 1)
    c_2 = hermite_coeff(n - 2)
    c = copy(c_1)
    if length(c_2) == length(c)
        push!(c, 0)
    end
    for i in 1:length(c_2)
        c[i + 1] = c[i + 1] - (n - 1) * c_2[i] / 2
    end
    @assert length(_hermit_coeff_cache) == n
    push!(_hermit_coeff_cache, c)
    return c
end

# Integral of x^2n exp(-x^2 / a)
function gaussian_int(n, a::T) where T
    return √a * a^n * gamma(T(n) + 0.5)
end

function poly_multiply(c1::AbstractVector{T1}, c2::AbstractVector{T2}) where {T1, T2}
    n1 = length(c1)
    n2 = length(c2)
    c = zeros(promote_type(T1, T2), n1 + n2 - 1)
    @inbounds for i in 1:n1
        v1 = c1[i]
        @simd for j in 1:n2
            c[i + j - 1] += v1 * c2[j]
        end
    end
    return c
end

# calculate c[end - i] * x^i
function poly_apply(c::AbstractVector{T1}, x::T2) where {T1,T2}
    n = length(c)
    T = promote_type(T1, T2)
    res = uninit_ary(Vector{T}, n)
    v::T = one(T)
    for i in 0:(n - 1)
        res[n - i] = c[n - i] * v
        v = v * x
    end
    return res
end

function wavefunction_overlap(n1::Integer, n2::Integer, m1::Integer, m2::Integer,
                              z1::T1, z2::T2) where {T1,T2}
    T = float(promote_type(T1, T2))
    if (n1 + n2 + m1 + m2) % 2 != 0
        return T(0)
    end
    c0 = 1 / π / z1 / z2
    # log of √(2^(n1 + n2 + m1 + m2) / (n1!n2!m1!m2!))
    lc1 = log(T(2)) / 2 * (n1 + n2 + m1 + m2)
    lc2 = (lgamma(T(n1 + 1)) + lgamma(T(n2 + 1)) + lgamma(T(m1 + 1)) + lgamma(T(m2 + 1))) / 2
    # This is the prefactor in front of the integral.
    c1 = c0 * exp(lc1 - lc2)
    z1² = z1^2
    z2² = z2^2
    # The a in the exponent of the Gaussian integral
    # 1 / a = 1 / z1^2 + 1 / z2^2 = (z1^2 + z2^2) / (z1^2 * z2^2)
    # a = (z1^2 * z2^2) / (z1^2 + z2^2)
    ga = (z1² * z2²) / (z1² + z2²)
    # The coefficients for the products of two hermite polynomial
    hc1 = poly_multiply(hermite_coeff(n1), hermite_coeff(n2)) # hc1[1] is n1+n2 power
    hc2 = poly_multiply(hermite_coeff(m1), hermite_coeff(m2)) # hc2[1] is m1+m2 power
    # The power of the last term
    p0_1 = (n1 + n2) - 2 * (length(hc1) - 1)
    p0_2 = (m1 + m2) - 2 * (length(hc2) - 1)
    c1 = c1 / z1^p0_1 / z2^p0_2

    # hc3[1] is n1 + n2 + m1 + m2 power
    hc3 = poly_multiply(poly_apply(hc1, 1 / z1²), poly_apply(hc2, 1 / z2²))

    n0 = (n1 + n2 + m1 + m2) ÷ 2

    intv = T(0.0)
    for i in 1:length(hc3)
        intv += gaussian_int(n0 + 1 - i, ga) * hc3[i]
    end
    return c1 * intv
end
