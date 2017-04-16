#!/usr/bin/julia

# Compute Rabi flopping with the present of decay terms
# The Hamiltonian is assumed to be time independent and the Rabi drive is on-resonance

using NaCsCalc.Utils: thread_rng, trand

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

struct RabiDecayParams{T}
    Γ₁::T
    Γ₂::T
    Δ::T
    Ω::T
    Γ::T
    Δ²::T
    Ω²::T
    Γ²::T
    Ω′²::T # Or Δ′²
    Ω′::T # Or Δ′
    ΔΩ′::T # Or ΔΔ′
    c1::T
    c2::T
    c3::T
    overdamp::Bool
    @inline function RabiDecayParams{T}(Ω::T, Γ₁::T, Γ₂::T) where T
        Δ = (Γ₁ - Γ₂) / 2
        Γ = (Γ₁ + Γ₂) / 2
        Δ² = Δ^2
        Ω² = Ω^2
        Γ² = Γ^2
        overdamp = Ω² < Δ²
        if overdamp
            Ω′² = Δ² - Ω²
        else
            Ω′² = Ω² - Δ²
        end
        Ω′ = sqrt(Ω′²)
        ΔΩ′ = Δ * Ω′
        if overdamp
            c1 = 2 * Γ * Ω²
            c2 = Δ * (Γ - Ω′) * (Ω′ - Δ)
            c3 = -Δ * (Γ + Ω′) * (Δ + Ω′)
        else
            c1 = Γ * Ω²
            c2 = ΔΩ′ * Γ₁
            c3 = Δ * muladd(Γ, -Δ, Ω′²)
        end
        return new(Γ₁, Γ₂, Δ, Ω, Γ, Δ², Ω², Γ², Ω′², Ω′, ΔΩ′, c1, c2, c3, overdamp)
    end
end

function propagate_2states_underdamp(params::RabiDecayParams{T}, tmax, rd) where T
    r = T(rand(rd)) * params.Ω′²
    # Now find the t for which `ψ^2(t) * Ω′² = r`.
    # First check if `ψ^2(tmax) * Ω′² > r`
    t::T = tmax
    # Tolerance
    yδ = max(T(2e-7), eps(T) * 10) * params.Ω′²
    if params.Ω′² - r <= yδ
        return zero(T), 1, (one(T), zero(T))
    end
    Ω′t::T = 0
    sinΩ′t::T = 0
    cosΩ′t::T = 0
    expΓt::T = 0
    ψ²::T = 0
    # ψ² * Ω′² is bound between exp(-Γt) * Ω * (Ω ± Δ)
    # Use this to compute a better bounds
    absΔ = abs(params.Δ)
    exp_lo = r / params.Ω / (params.Ω + absΔ)
    @fastmath thi::T = -log(exp_lo) / params.Γ
    if !(thi < t)
        thi = t
        Ω′t = params.Ω′ * t
        sinΩ′t, cosΩ′t = sincos(Ω′t)
        @fastmath expΓt = exp(-params.Γ * t)
        ψ² = -expΓt * muladd(params.ΔΩ′, sinΩ′t, muladd(params.Δ², cosΩ′t, -params.Ω²))
        if ψ² > r
            # No decay happened, return the wave function at tmax
            Ω′t_2 = Ω′t / 2
            sinΩ′t_2, cosΩ′t_2 = sincos(Ω′t_2)
            ψ1 = muladd(params.Ω′, cosΩ′t_2, -params.Δ * sinΩ′t_2)
            ψ2 = params.Ω * sinΩ′t_2
            @fastmath factor = sqrt(expΓt / ψ²)
            return t, 0, (ψ1 * factor, ψ2 * factor)
        elseif r - ψ² <= yδ
            @goto ret
        end
    end
    exp_hi = r / params.Ω / (params.Ω - absΔ)
    @fastmath tlo::T = -log(exp_hi) / params.Γ
    if !(tlo > 0)
        tlo = 0
    end
    # Find the root using a combination of newton's method and bisecting.
    # The bisection is to avoid newton's method overshotting since the function we want
    # to solve is known to have partial oscillation.
    tthresh = min(T(2 / params.Ω′), T(tmax))
    first_loop = true
    while thi - tlo > 5 * eps(T) * thi || first_loop
        if thi - tlo < tthresh && !first_loop
            # Try Newton's method
            diff = -expΓt * muladd(params.c3, cosΩ′t, -muladd(params.c2, sinΩ′t, -params.c1))
            δ1 = (ψ² - r)
            t2 = t - δ1 / diff
            if tlo < t2 < thi
                δt = thi - tlo
                t = t2
                Ω′t = params.Ω′ * t
                sinΩ′t, cosΩ′t = sincos(Ω′t)
                @fastmath expΓt = exp(-params.Γ * t)
                ψ² = -expΓt * muladd(params.ΔΩ′, sinΩ′t, muladd(params.Δ², cosΩ′t, -params.Ω²))
                δ2 = ψ² - r
                # We've already computed this point, even if it's not enough this time,
                # might as well use it to narrow down the range
                if δ2 < 0
                    δ2 = -δ2
                    δ2 < yδ && @goto ret
                    thi = t
                else
                    δ2 < yδ && @goto ret
                    tlo = t
                end
                δt2 = thi - tlo
                # If we shrink either x or y by at least √2, keep using the Newton's method
                # otherwise, shrink it further with bisect.
                if δt2 * 1.4 <= δt || δ2 * 1.4 < abs(δ1)
                    continue
                end
            end
        end
        first_loop = false
        t = (thi + tlo) / 2
        Ω′t = params.Ω′ * t
        sinΩ′t, cosΩ′t = sincos(Ω′t)
        @fastmath expΓt = exp(-params.Γ * t)
        ψ² = -expΓt * muladd(params.ΔΩ′, sinΩ′t, muladd(params.Δ², cosΩ′t, -params.Ω²))
        δ2 = ψ² - r
        if δ2 < 0
            -δ2 < yδ && @goto ret
            thi = t
        else
            δ2 < yδ && @goto ret
            tlo = t
        end
    end

    @label ret
    # Whether we decay through state 1 or 2 here depends on the instantaneous decay rate
    # at `t`
    rtotal = 2 * muladd(params.c3, cosΩ′t, -muladd(params.c2, sinΩ′t, -params.c1))
    r2 = params.Γ₂ * params.Ω² * (1 - cosΩ′t)
    return t, rtotal * rand(rd) < r2 ? 2 : 1, (one(T), zero(T))
end

function propagate_2states_overdamp(params::RabiDecayParams{T}, tmax, rd) where T
    r0 = T(rand(rd))
    # Now find the t for which `ψ^2(t) = r0`.
    t::T = tmax
    # Tolerance
    y0δ = max(T(2e-7), eps(T) * 10)
    if 1 - r0 <= y0δ
        return zero(T), 1, (one(T), zero(T))
    end
    r = r0 * 2 * params.Ω′²
    yδ = y0δ * 2 * params.Ω′²
    ΔmΔ′ = params.Δ - params.Ω′
    ΔpΔ′ = params.Δ + params.Ω′

    Δ′t::T = 0
    expΔ′t::T = 0
    exp_Δ′t::T = 0
    expΓt::T = 0
    ψ²::T = 0
    # ψ² is bound between exp(-(Γ ± Δ)t)
    absΔ = abs(params.Δ)
    thi::T = -@fastmath(log(r0)) / (params.Γ - absΔ)
    if !(thi < t)
        thi = t
        Δ′t = params.Ω′ * t
        @fastmath expΔ′t = exp(Δ′t)
        exp_Δ′t = 1 / expΔ′t
        @fastmath expΓt = exp(-params.Γ * t)
        ψ² = expΓt * muladd(params.Δ, muladd(ΔmΔ′, expΔ′t, ΔpΔ′ * exp_Δ′t), -2 * params.Ω²)
        if ψ² > r
            # No decay happened, return the wave function at tmax
            Δ′t_2 = Δ′t / 2
            @fastmath expΔ′t_2 = exp(Δ′t_2)
            exp_Δ′t_2 = 1 / expΔ′t_2
            ψ1 = muladd(ΔpΔ′, exp_Δ′t_2, -ΔmΔ′ * expΔ′t_2)
            ψ2 = params.Ω * (expΔ′t_2 - exp_Δ′t_2)
            @fastmath factor = sqrt(expΓt / ψ² / 2)
            return t, 0, (ψ1 * factor, ψ2 * factor)
        elseif r - ψ² <= yδ
            @goto ret
        end
    end
    tlo::T = -@fastmath(log(r0)) / (params.Γ + absΔ)
    t = tlo
    Δ′t = params.Ω′ * t
    @fastmath expΔ′t = exp(Δ′t)
    exp_Δ′t = 1 / expΔ′t
    @fastmath expΓt = exp(-params.Γ * t)
    ψ² = expΓt * muladd(params.Δ, muladd(ΔmΔ′, expΔ′t, ΔpΔ′ * exp_Δ′t), -2 * params.Ω²)
    # Find the root using a combination of newton's method and bisecting.
    # The bisection is to avoid newton's method overshotting.
    # It shouldn't be needed in this case since the function does not oscillate but is
    # included just as a fallback.
    while thi - tlo > 5 * eps(T) * thi
        # Try Newton's method
        diff = expΓt * muladd(params.c3, exp_Δ′t, muladd(params.c2, expΔ′t, params.c1))
        δ1 = (ψ² - r)
        t2 = t - δ1 / diff
        if tlo < t2 < thi
            δt = thi - tlo
            t = t2
            Δ′t = params.Ω′ * t
            @fastmath expΔ′t = exp(Δ′t)
            exp_Δ′t = 1 / expΔ′t
            @fastmath expΓt = exp(-params.Γ * t)
            ψ² = expΓt * muladd(params.Δ, muladd(ΔmΔ′, expΔ′t, ΔpΔ′ * exp_Δ′t), -2 * params.Ω²)
            δ2 = ψ² - r
            # We've already computed this point, even if it's not enough this time,
            # might as well use it to narrow down the range
            if δ2 < 0
                δ2 = -δ2
                δ2 < yδ && @goto ret
                thi = t
            else
                δ2 < yδ && @goto ret
                tlo = t
            end
            δt2 = thi - tlo
            # If we shrink either x or y by at least √2, keep using the Newton's method
            # otherwise, shrink it further with bisect.
            if δt2 * 1.4 <= δt || δ2 * 1.4 < abs(δ1)
                continue
            end
        end
        t = (thi + tlo) / 2
        Δ′t = params.Ω′ * t
        @fastmath expΔ′t = exp(Δ′t)
        exp_Δ′t = 1 / expΔ′t
        @fastmath expΓt = exp(-params.Γ * t)
        ψ² = expΓt * muladd(params.Δ, muladd(ΔmΔ′, expΔ′t, ΔpΔ′ * exp_Δ′t), -2 * params.Ω²)
        δ2 = ψ² - r
        if δ2 < 0
            -δ2 < yδ && @goto ret
            thi = t
        else
            δ2 < yδ && @goto ret
            tlo = t
        end
    end

    @label ret
    # Whether we decay through state 1 or 2 here depends on the instantaneous decay rate
    # at `t`
    rtotal = -2 * muladd(params.c3, exp_Δ′t, muladd(params.c2, expΔ′t, params.c1))
    r2 = params.Γ₂ * params.Ω² * (exp_Δ′t + expΔ′t - 2)
    return t, rtotal * rand(rd) < r2 ? 2 : 1, (one(T), zero(T))
end

"""
    propagate_2states(params::RabiDecayParams, tmax, rd) -> (t, i, ψ)

The atom start in state 1 and is doing Rabi flopping with (angular) Rabi frequency `Ω`
to state 2. The state 1(2) has a decay rate of `Γ₁`(`Γ₂`). Propagate this system under the
modified Hamiltonian for quantum jump method for at most `tmax` or until a decay happens.

Returns the actual time propagated `t`.
If a decay has happend, `t` is the decay time. `i` (either `1` or `2`) is the state from which
the decay occurs. `ψ` is unused.
If no decay happens, `t == tmax`, `i == 0`, `ψ` is a tuple of the wavefunctions
"""
@inline propagate_2states(params::RabiDecayParams, tmax, rd) = if params.overdamp
    propagate_2states_overdamp(params, tmax, rd)
else
    propagate_2states_underdamp(params, tmax, rd)
end

"""
    propagate(Ω, Γ, rates, tmax, rd=thread_rng()) -> ψ

The atom start in state 1 and is doing Rabi flopping with (angular) Rabi frequency `Ω`
to state 2.
The decay rates and matrix are given by `Γ` and `rates` which are assumed to satisfy
`ratesⱼ = ∑ᵢΓᵢⱼ`.
Propagate this system using the quantum jump method once.

Returns the wave function after the propagation.
"""
function propagate(Ω::T, Γ::AbstractMatrix{T}, rates::AbstractVector{T},
                   _tmax, rd=thread_rng()) where T
    Γ₁, Γ₂ = rates
    tmax::T = _tmax
    params1 = RabiDecayParams{T}(Ω, Γ₁, Γ₂)
    params2 = RabiDecayParams{T}(Ω, Γ₂, Γ₁)
    flipped = false
    @inbounds while tmax > 0
        t, idx, ψ = propagate_2states(params1, tmax, rd)
        tmax -= t
        if idx == 0
            if flipped
                return ψ[2], ψ[1]
            else
                return ψ
            end
        end
        if flipped
            idx = 3 - idx
        end
        if rand(rd) * rates[idx] > Γ[1, idx]
            if !flipped
                flipped = true
                params1, params2 = params2, params1
            end
        else
            if flipped
                flipped = false
                params1, params2 = params2, params1
            end
        end
    end
    return flipped ? (zero(T), one(T)) : (one(T), zero(T))
end

"""
    average(Ω, Γ, rates, tmax, n, rd=thread_rng()) -> ψ, σψ

The atom start in state 1 and is doing Rabi flopping with (angular) Rabi frequency `Ω`
to state 2.
The decay rates and matrix are given by `Γ` and `rates` which are assumed to satisfy
`ratesⱼ = ∑ᵢΓᵢⱼ`.
Propagate this system using the quantum jump method by `n` times.

Returns the averaged probability distribution and its uncertainty after the propagation.
"""
function average(Ω::T, Γ::AbstractMatrix{T}, rates::AbstractVector{T}, tmax, n,
                 rd=thread_rng()) where T
    # Always use Float64 for the sum so that the order of summing does not matter as much
    # we could also be fancier and use the real type with the optimum summing order but
    # using Float64 for sum is cheap and the easiest way to implement.
    sumϕ1 = 0.0
    sumϕ2 = 0.0
    sumϕ²1 = 0.0
    sumϕ²2 = 0.0
    @inbounds for i in 1:n
        ψ = propagate(Ω, Γ, rates, tmax, rd)
        ψ1 = abs2(ψ[1])
        ψ2 = abs2(ψ[2])
        sumϕ1 += ψ1
        sumϕ²1 += ψ1^2
        sumϕ2 += ψ2
        sumϕ²2 += ψ2^2
    end
    sqrtn = sqrt(n - 1)
    avgϕ1 = sumϕ1 / n
    avgϕ²1 = sumϕ²1 / n
    σϕ1 = (avgϕ²1 - avgϕ1^2) / sqrtn
    avgϕ2 = sumϕ2 / n
    avgϕ²2 = sumϕ²2 / n
    σϕ2 = (avgϕ²2 - avgϕ2^2) / sqrtn
    return (avgϕ1, avgϕ2), (σϕ1, σϕ2)
end

function Γ_to_rates(Γ)
    len = size(Γ, 1)
    rates = Vector{eltype(Γ)}(len)
    for i in 1:len
        local s = zero(eltype(Γ))
        for j in 1:len
            s += Γ[j, i]
        end
        rates[i] = s
    end
    return rates
end

ϕ = [1.0, 0.0]
i1 = 1
i2 = 2
const δt = 1e-8

using PyPlot
pts = 0:200:10000
res = Vector{Float64}(length(pts))
unc = Vector{Float64}(length(pts))

function f(pts, Γ, ϕ, i1, i2, Ω, res, unc, color)
    len = length(ϕ)
    rates = Γ_to_rates(Γ)
    Γ₁, Γ₂ = rates
    params1 = RabiDecayParams{Float64}(Ω, Γ₁, Γ₂)
    dump(params1)
    # res .= 0
    # unc .= 0
    # @time Threads.@threads for i in 1:length(pts)
    #     local a, s
    #     a, s = average(Ω, Γ, rates, δt * pts[i], 10000)
    #     res[i] = a[1]
    #     unc[i] = s[1]
    # end
    # errorbar(pts * δt, res, unc, fmt="-", label="0", color=color)
    Ω32 = Float32(Ω)
    Γ32 = Float32.(Γ)
    rates32 = Float32.(rates)
    δt32 = Float32(δt)
    res .= 0
    unc .= 0
    @time Threads.@threads for i in 1:length(pts)
        local a, s
        a, s = average(Ω32, Γ32, rates32, δt32 * pts[i], 10000)
        res[i] = a[1]
        unc[i] = s[1]
    end
    errorbar(pts * δt, res, unc, fmt="^-", label="0", color=color)
end
Ω = 2π * 20e3
Γ = [2e4 0
      0 3e4] * 4
f(pts, Γ, ϕ, i1, i2, Ω, res, unc, "blue")
Ω = 2π * 20e3
Γ = [0 1e4
      3e4 0] * 4
f(pts, Γ, ϕ, i1, i2, Ω, res, unc, "red")
Ω = 2π * 2e3
Γ = [2e4 0
      0 3e4] * 4
f(pts, Γ, ϕ, i1, i2, Ω, res, unc, "green")
Ω = 2π * 2e3
Γ = [0 1e4
      3e4 0] * 4
f(pts, Γ, ϕ, i1, i2, Ω, res, unc, "orange")
# Γ = [0 3e4
#       3e4 0] * 2
# f(pts, Γ, ϕ, i1, i2, Ω, res, unc, "red")
# ts = δt * pts
# function y(t, _Ω)
#     _Γ = 3e4 * 2
#     Ω = √(_Ω^2 + 2 * _Γ^2)
#     Γ = 3_Γ
#     Ω′ = √(Ω^2 - Γ^2 / 4)
#     -exp(-Γ * t / 2) * (Γ / 2 / Ω′ * sin(Ω′ * t) + cos(Ω′ * t))
# end
# plot(ts, (1 .- y.(ts, Ω)) ./ 2, color="orange")
# Γ = [3e4 0
#       0 3e4] * 2
# f(pts, Γ, ϕ, i1, i2, Ω, res, unc, "blue")
# ts = δt * pts
# function y(t, Ω)
#     Γ = 3e4 * 2
#     Ω′ = √(Ω^2 - Γ^2 / 4)
#     -exp(-Γ * t / 2) * (Γ / 2 / Ω′ * sin(Ω′ * t) + cos(Ω′ * t))
# end
# plot(ts, (1 .- y.(ts, Ω)) ./ 2, color="orange")
# Γ = [0 4e4
#       2e4 0] * 2
# f(pts, Γ, ϕ, i1, i2, Ω, res, unc, "cyan")
# Γ = [0 6e4
#       0e4 0] * 2
# f(pts, Γ, ϕ, i1, i2, Ω, res, unc, "orange")
# Γ = [0 3e4
#       3e4 0] * 2
# f(pts, Γ, ϕ, i1, i2, 0.1e3, res, unc, "blue")

legend()
grid()

# @show propagate(Ω, Γ, rates, 2e-9 * 1000, Base.Random.GLOBAL_RNG)

function f2(Ω, Γ)
    rates = Γ_to_rates(Γ)
    propagate(Ω, Γ, rates, 0.11e-3)
    @time Threads.@threads for i in 1:100000000
        propagate(Ω, Γ, rates, 0.11e-3)
    end
end
f2(Ω, Γ)

show()

# @code_native propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
