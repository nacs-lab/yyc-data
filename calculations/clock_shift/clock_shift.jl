#

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
    assert(length(_hermit_coeff_cache) == n)
    push!(_hermit_coeff_cache, c)
    return c
end

# Integral of x^2n exp(-x^2 / a)
function gaussian_int(n, a)
    return √a * a^n * gamma(n + 0.5)
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

# calculate c[end - i] * x^i * a
function poly_apply(c::AbstractVector{T1}, x::T2, a::T3=1) where {T1,T2,T3}
    n = length(c)
    T = promote_type(T1, T2, T3)
    res = Vector{T}(n)
    v::T = one(T)
    for i in 0:(n - 1)
        res[n - i] = c[n - i] * v * a
        v = v * x
    end
    return res
end
