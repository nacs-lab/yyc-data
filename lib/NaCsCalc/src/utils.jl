#!/usr/bin/julia

module Utils

# LGPLv3 implementation from libstdc++

function poly_laguerre_large_n{Tp}(n::Integer, α, x::Tp)
    a::Tp = -n
    b::Tp = α + 1
    η::Tp = 2b - 4a
    cos²th = x / η
    sin²th = 1 - cos²th
    costh = @fastmath sqrt(cos²th)
    th = @fastmath acos(costh)
    pre_h = (Tp(π / 2)^2) * η * η * cos²th * sin²th
    lg_b = @fastmath lgamma(Tp(n + b))
    lnfact = @fastmath lgamma(Tp(n + 1))
    pre_term1 = @fastmath Tp(0.5) * (1 - b) * log(Tp(0.25) * x * η)
    pre_term2 = @fastmath Tp(0.25) * log(pre_h)
    lnpre = lg_b - lnfact + Tp(0.5) * x + pre_term1 - pre_term2

    th2 = 2 * th
    sin2th = @fastmath 2 * costh * sqrt(sin²th)

    # From libstdc++
    ser_term1 = 0 # @fastmath sinpi(a)
    # This might be off by a minus sign or sth like that.
    # Evaluating at `10000001, 10, 1.2` gives the wrong result
    ser_term2 = @fastmath sin(Tp(0.25) * η * (th2 - sin2th) + Tp(π / 4))

    ser::Tp = ser_term1 + ser_term2
    return @fastmath exp(lnpre) * ser
end

function poly_laguerre_hyperg{Tp}(n::Integer, α, x::Tp)
    b::Tp = Tp(α) + 1
    mx = -x
    tc_sgn::Tp = x < 0 ? 1 : ((n % 2 == 1) ? -1 : 1)
    # Get |x|^n/n!
    tc::Tp = 1
    ax = abs(x)
    for k in 1:n
        tc *= ax / k
    end
    term::Tp = tc * tc_sgn
    _sum::Tp = term
    for k in (n - 1):-1:0
        term *= ((b + Tp(k)) / Tp(n - k)) * Tp(k + 1) / mx
        _sum += term
    end
    return _sum
end

function poly_laguerre_recursion{Tp}(n::Integer, α, x::Tp)
    # Compute l_0.
    l_0::Tp = 1
    n == 0 && return l_0

    # Compute l_1^alpha.
    l_1::Tp = -x + 1 + α
    n == 1 && return l_1

    # Compute l_n^alpha by recursion on n.
    l_n′::Tp = l_0
    l_n::Tp = l_1
    b::Tp = α - 1
    a::Tp = b - x
    @fastmath for nn in 2:n
        fnn = Tp(nn)
        l1 = muladd(a, l_n, -b * l_n′)
        l2 = muladd(2, l_n, -l_n′)
        l_n, l_n′ = l1 / fnn + l2, l_n
    end
    return l_n
end

function genlaguerre{Tp<:AbstractFloat}(n::Integer, α, x::Tp)::Tp
    if x < 0
        throw(DomainError())
    elseif isnan(α) || isnan(x)
        # Return NaN on NaN input.
        return NaN
    elseif n == 0
        return 1
    elseif n == 1
        return Tp(1) + Tp(α) - x
    elseif x == 0
        prod::Tp = α + 1
        for k in 2:n
            prod *= Tp(α + k) / Tp(k)
        end
        return prod
    elseif n > 10000000 && α > -1 && x < 2 * (α + 1) + 4n
        return poly_laguerre_large_n(n, α, x)
    elseif α >= 0 || (x > 0 && α < -(n + 1))
        return poly_laguerre_recursion(n, α, x)
    else
        return poly_laguerre_hyperg(n, α, x)
    end
end
genlaguerre(n::Integer, α, x) = genlaguerre(n, α, float(x))

function binomial_estimate{T<:AbstractFloat}(x, n, z::T=1.0)
    if n <= 0
        return T(0.5), T(0.5)
    end
    p = T(x / n)
    z² = z^2
    z²n = z² / n
    p′::T = (p + z²n / 2) / (1 + z²n)
    unc::T = sqrt(p * (1 - p) / n + z² / 4 / n^2) / (1 + z²n)
    return p′, unc
end

@inline function sincos(v::Float64)
    Base.llvmcall("""
    %f = bitcast i8 *%1 to void (double, double *, double *)*
    %pres = alloca [2 x double]
    %p1 = getelementptr inbounds [2 x double], [2 x double]* %pres, i64 0, i64 0
    %p2 = getelementptr inbounds [2 x double], [2 x double]* %pres, i64 0, i64 1
    call void %f(double %0, double *nocapture noalias %p1, double *nocapture noalias %p2)
    %res = load [2 x double], [2 x double]* %pres
    ret [2 x double] %res
    """, Tuple{Float64,Float64}, Tuple{Float64,Ptr{Void}}, v, cglobal((:sincos, Base.libm_name)))
end

@inline function sincos(v::Float32)
    Base.llvmcall("""
    %f = bitcast i8 *%1 to void (float, float *, float *)*
    %pres = alloca [2 x float]
    %p1 = getelementptr inbounds [2 x float], [2 x float]* %pres, i64 0, i64 0
    %p2 = getelementptr inbounds [2 x float], [2 x float]* %pres, i64 0, i64 1
    call void %f(float %0, float *nocapture noalias %p1, float *nocapture noalias %p2)
    %res = load [2 x float], [2 x float]* %pres
    ret [2 x float] %res
    """, Tuple{Float32,Float32}, Tuple{Float32,Ptr{Void}}, v, cglobal((:sincosf, Base.libm_name)))
end

@inline function sincos(v)
    @fastmath (sin(v), cos(v))
end

const ThreadRNG = MersenneTwister[]
function __init__()
    # Allocate the random number generator on the thread's own heap
    # instead of the master thread heap to minimize memory conflict
    nth = Threads.nthreads()
    resize!(ThreadRNG, nth)
    init_rng = function ()
        tid = Threads.threadid()
        ThreadRNG[tid] = MersenneTwister(0)
    end
    ccall(:jl_threading_run, Ref{Void}, (Any,), init_rng)
end

@inline function thread_rng()
    # Bypass bounds check and NULL check
    Base.llvmcall(
        """
        %p = load i8**, i8*** %0
        ret i8** %p
        """, MersenneTwister, Tuple{Ptr{Ptr{Ptr{Void}}}},
        Ptr{Ptr{Ptr{Void}}}(pointer(ThreadRNG) + (Threads.threadid() - 1) * sizeof(Int)))
end

@inline trand() = rand(thread_rng())

end
