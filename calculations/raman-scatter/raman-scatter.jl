#!/usr/bin/julia

# Compute Rabi flopping with the present of decay terms
# The Hamiltonian is assumed to be time independent and the Rabi drive is on-resonance

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

@inline function sincosh(v)
    @fastmath expv = exp(v)
    exp_v = 1 / expv
    sinhv = (expv - exp_v) / 2
    coshv = (expv + exp_v) / 2
    return sinhv, coshv
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
        c1 = Γ * Ω²
        c2 = ΔΩ′ * Γ₁
        c3 = Δ * muladd(Γ, -Δ, Ω′²)
        return new(Γ₁, Γ₂, Δ, Ω, Γ, Δ², Ω², Γ², Ω′², Ω′, ΔΩ′, c1, c2, c3, overdamp)
    end
end

"""
    propagate_2states(params::RabiDecayParams, tmax) -> (t, i, ψ)

The atom start in state 1 and is doing Rabi flopping with (angular) Rabi frequency `Ω`
to state 2. The state 1(2) has a decay rate of `Γ₁`(`Γ₂`). Propagate this system under the
modified Hamiltonian for quantum jump method for at most `tmax` or until a decay happens.

Returns the actual time propagated `t`.
If a decay has happend, `t` is the decay time. `i` (either `1` or `2`) is the state from which
the decay occurs. `ψ` is unused.
If no decay happens, `t == tmax`, `i == 0`, `ψ` is a tuple of the wavefunctions
"""
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
    r = T(rand(rd)) * params.Ω′²
    # Now find the t for which `ψ^2(t) * Δ′² = r`.
    # First check if `ψ^2(tmax) * Δ′² > r`
    t::T = tmax
    # Tolerance
    yδ = max(T(2e-7), eps(T) * 10) * params.Ω′²
    if params.Ω′² - r <= yδ
        return zero(T), 1, (one(T), zero(T))
    end
    Δ′t::T = params.Ω′ * t
    sinhΔ′t::T = 0
    coshΔ′t::T = 0
    sinhΔ′t, coshΔ′t = sincosh(Δ′t)
    @fastmath expΓt::T = exp(-params.Γ * t)
    ψ²::T = expΓt * muladd(params.Δ², coshΔ′t, -muladd(params.ΔΩ′, sinhΔ′t, params.Ω²))
    if ψ² > r
        # No decay happened, return the wave function at tmax
        Δ′t_2 = Δ′t / 2
        sinhΔ′t_2, coshΔ′t_2 = sincosh(Δ′t_2)
        ψ1 = muladd(params.Ω′, coshΔ′t_2, -params.Δ * sinhΔ′t_2)
        ψ2 = params.Ω * sinhΔ′t_2
        @fastmath factor = sqrt(expΓt / ψ²)
        return t, 0, (ψ1 * factor, ψ2 * factor)
    elseif r - ψ² <= yδ
        @goto ret
    end
    thi::T = t
    tlo::T = 0
    if r^2 > ψ² * params.Ω′²
        # It's closer to t=0, start from there
        t = 0
        Δ′t = 0
        sinhΔ′t = 0
        coshΔ′t = 1
        expΓt = 1
        ψ² = params.Ω′²
    end
    # Find the root using a combination of newton's method and bisecting.
    # The bisection is to avoid newton's method overshotting.
    # It shouldn't be needed in this case since the function does not oscillate but is
    # included just as a fallback.
    while thi - tlo > 5 * eps(T) * thi
        # Try Newton's method
        diff = expΓt * muladd(params.c3, coshΔ′t, muladd(params.c2, sinhΔ′t, params.c1))
        δ1 = (ψ² - r)
        t2 = t - δ1 / diff
        if tlo < t2 < thi
            δt = thi - tlo
            t = t2
            Δ′t = params.Ω′ * t
            sinhΔ′t, coshΔ′t = sincosh(Δ′t)
            @fastmath expΓt = exp(-params.Γ * t)
            ψ² = expΓt * muladd(params.Δ², coshΔ′t, -muladd(params.ΔΩ′, sinhΔ′t, params.Ω²))
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
        sinhΔ′t, coshΔ′t = sincosh(Δ′t)
        @fastmath expΓt = exp(-params.Γ * t)
        ψ² = expΓt * muladd(params.Δ², coshΔ′t, -muladd(params.ΔΩ′, sinhΔ′t, params.Ω²))
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
    rtotal = -2 * muladd(params.c3, coshΔ′t, muladd(params.c2, sinhΔ′t, params.c1))
    r2 = params.Γ₂ * params.Ω² * (coshΔ′t - 1)
    return t, rtotal * rand(rd) < r2 ? 2 : 1, (one(T), zero(T))
end

function propagate2(Ω::T, Γ::AbstractMatrix{T}, rates::AbstractVector{T}, _tmax, rd) where T
    Γ₁, Γ₂ = rates
    tmax::T = _tmax
    params1 = RabiDecayParams{T}(Ω, Γ₁, Γ₂)
    params2 = RabiDecayParams{T}(Ω, Γ₂, Γ₁)
    overdamp = params1.overdamp
    flipped = false
    @inbounds while tmax > 0
        t, idx, ψ = if overdamp
            propagate_2states_overdamp(params1, tmax, rd)
        else
            propagate_2states_underdamp(params1, tmax, rd)
        end
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

function average2(Ω::T, Γ::AbstractMatrix{T}, rates::AbstractVector{T}, tmax, n, rd) where T
    sumϕ1 = zero(T)
    sumϕ2 = zero(T)
    sumϕ²1 = zero(T)
    sumϕ²2 = zero(T)
    @inbounds for i in 1:n
        ψ = propagate2(Ω, Γ, rates, tmax, rd)
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

# Propagate the state `ϕ` by `nstep` number of time steps each of length `δt`
# The decay matrix is `Γ`, the element `(i,j)` of this matrix represents the decay rate
# from state `j` to state `i`.
# The Rabi flopping happens on resonant from state `i1` to state `i2` with Rabi frequency `Ω`
function propagate(ϕ, δt, nstep, Γ, i1, i2, Ω, nstate, rates, buff1, rd)
    # Do propagation of coherent and dissipative part separately.
    # Use precise exponentiation for the (potentially fast) Rabi flopping and use linear
    # approximation for the decay probability.
    # Compute the decay probability in the middle of the time step to improve precision.
    if nstep == 0
        return ϕ
    end

    # Half step
    @fastmath sinΩt_2 = sin(Ω * δt / 4)
    @fastmath cosΩt_2 = cos(Ω * δt / 4)
    # Full step
    @fastmath sinΩt = sin(Ω * δt / 2)
    @fastmath cosΩt = cos(Ω * δt / 2)
    # Decay rate
    T = real(eltype(ϕ))
    s::T = 0
    for i in 1:nstate
        s += abs2(ϕ[i])
    end

    # First compute the shifted wavefunction. We'll propagate full steps from here
    # until the last step.
    # This also checks bound so we can skip bound checking later
    ϕ1 = ϕ[i1]
    ϕ2 = ϕ[i2]
    ϕ1, ϕ2 = ϕ1 * cosΩt_2 + ϕ2 * sinΩt_2, ϕ2 * cosΩt_2 - ϕ1 * sinΩt_2
    ϕ[i1] = ϕ1
    ϕ[i2] = ϕ2
    @inbounds for i in 1:(nstep - 1)
        # We are at `t = (i - 1/2) * δt`
        # Decide if we need to decay first
        @simd for j in 1:nstate
            buff1[j] = abs2(ϕ[j]) * rates[j] * δt
        end
        # buff1 .= abs2.(ϕ) .* rates .* δt
        r = rand(rd) * s
        decay = 0
        for j in 1:nstate
            r -= buff1[j]
            if r < 0
                decay = j
                break
            end
        end
        if decay != 0
            # decay happend from state `decay`
            # Decide which state it'll fall into
            r = rand(rd) * rates[decay]
            local k = 0
            for k in 1:nstate
                r -= Γ[k, decay]
                if r < 0
                    break
                end
            end
            s = 1
            ϕ .= 0
            if i1 == k
                ϕ[i1] = cosΩt_2
                ϕ[i2] = -sinΩt_2
            elseif i2 == k
                ϕ[i1] = sinΩt_2
                ϕ[i2] = cosΩt_2
            else
                ϕ[k] = 1
            end
            continue
        end
        # No decay, cool. Propagate by a full time step to `t = (i + 1/2) * δt`
        # First propagate the decay term
        @fastmath factor = 1 / 2 / sqrt(s)
        s = 0 # normalization
        @simd for j in 1:nstate
            v = ϕ[j] * (1 - rates[j] * δt / 2) * factor
            s += abs2(v)
            ϕ[j] = v
        end
        # Then propagate the drive
        ϕ1 = ϕ[i1]
        ϕ2 = ϕ[i2]
        ϕ1, ϕ2 = ϕ1 * cosΩt + ϕ2 * sinΩt, ϕ2 * cosΩt - ϕ1 * sinΩt
        ϕ[i1] = ϕ1
        ϕ[i2] = ϕ2
    end
    ϕ1 = ϕ[i1]
    ϕ2 = ϕ[i2]
    ϕ1, ϕ2 = ϕ1 * cosΩt_2 + ϕ2 * sinΩt_2, ϕ2 * cosΩt_2 - ϕ1 * sinΩt_2
    ϕ[i1] = ϕ1
    ϕ[i2] = ϕ2
    ϕ ./= √(s)
    return ϕ
end

function average(ϕ0, δt, nstep, Γ, i1, i2, Ω, n, rates, rd)
    len = length(ϕ0)
    ϕ = similar(ϕ0)
    sumϕ = zeros(real(eltype(ϕ0)), len)
    sumϕ² = zeros(real(eltype(ϕ0)), len)
    T = real(eltype(ϕ))
    buff1 = Vector{T}(len)
    @inbounds for i in 1:n
        ϕ .= ϕ0
        propagate(ϕ, δt, nstep, Γ, i1, i2, Ω, len, rates, buff1, rd)
        for j in 1:len
            ϕi = abs2(ϕ[j])
            sumϕ[j] += ϕi
            sumϕ²[j] += ϕi^2
        end
    end
    @inbounds for i in 1:len
        avgϕi = sumϕ[i] / n
        avgϕ²i = sumϕ²[i] / n
        sumϕ[i] = avgϕi
        sumϕ²[i] = (avgϕ²i - avgϕi^2) / sqrt(n - 1)
    end
    return sumϕ, sumϕ²
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
Ω = 2π * 2e3
const δt = 1e-8

using PyPlot
pts = 0:20:6000
res = Vector{Float64}(length(pts))
unc = Vector{Float64}(length(pts))

function f(pts, Γ, ϕ, i1, i2, Ω, res, unc, color)
    len = length(ϕ)
    rates = Γ_to_rates(Γ)
    Γ₁, Γ₂ = rates
    params1 = RabiDecayParams{Float64}(Ω, Γ₁, Γ₂)
    dump(params1)
    rds = [MersenneTwister(0) for i in 1:Threads.nthreads()]
    res .= 0
    unc .= 0
    @time Threads.@threads for i in 1:length(pts)
        local a, s
        a, s = average(ϕ, δt, pts[i], Γ, i1, i2, Ω, 2000, rates,
                       rds[Threads.threadid()])
        res[i] = a[1]
        unc[i] = s[1]
    end
    errorbar(pts * δt, res, unc, fmt=".", label="1.1", color=color)
    # @time Threads.@threads for i in 1:length(pts)
    #     local a, s
    #     a, s = average(ϕ, δt / 100, pts[i] * 100, Γ, i1, i2, Ω, 10000,
    #                    rates, rds[Threads.threadid()])
    #     res[i] = a[1]
    #     unc[i] = s[1]
    # end
    # errorbar(pts * δt, res, unc, fmt=".", label="0.55")
    res .= 0
    unc .= 0
    @time Threads.@threads for i in 1:length(pts)
        local a, s
        a, s = average2(Ω, Γ, rates, δt * pts[i], 10000, rds[Threads.threadid()])
        res[i] = a[1]
        unc[i] = s[1]
    end
    errorbar(pts * δt, res, unc, fmt="-", label="0", color=color)
    Ω32 = Float32(Ω)
    Γ32 = Float32.(Γ)
    rates32 = Float32.(rates)
    δt32 = Float32(δt)
    res .= 0
    unc .= 0
    @time Threads.@threads for i in 1:length(pts)
        local a, s
        a, s = average2(Ω32, Γ32, rates32, δt32 * pts[i], 10000, rds[Threads.threadid()])
        res[i] = a[1]
        unc[i] = s[1]
    end
    errorbar(pts * δt, res, unc, fmt="^-", label="0", color=color)
end
Γ = [2e4 0
      0 3e4] * 4
f(pts, Γ, ϕ, i1, i2, Ω, res, unc, "blue")
Γ = [0 3e4
      2e4 0] * 4
f(pts, Γ, ϕ, i1, i2, Ω, res, unc, "red")
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
# for i in 1:npts
#     a, s = average(ϕ, δt, 10 * (i - 1), Γ, i1, i2, 0, 10000)
#     res[i] = a[1]
#     unc[i] = s[1]
# end
# errorbar(1:npts, res, unc)

legend()
grid()
show()

# @show propagate2(Ω, Γ, rates, 2e-9 * 1000, Base.Random.GLOBAL_RNG)
# # @time @show average(ϕ, 5.5e-7, 200, Γ, i1, i2, 0, 1)
# # @time @show average(ϕ, 2.75e-7, 40_000, Γ, i1, i2, Ω, 100)
# # @time @show average(ϕ, 1.1e-7, 100_000, Γ, i1, i2, Ω, 100)

function f2(Ω, Γ)
    rates = Γ_to_rates(Γ)
    rds = [MersenneTwister(0) for i in 1:Threads.nthreads()]
    propagate2(Ω, Γ, rates, 0.11e-3, rds[Threads.threadid()])
    @time Threads.@threads for i in 1:100000000
        propagate2(Ω, Γ, rates, 0.11e-3, rds[Threads.threadid()])
    end
end
f2(Ω, Γ)

# @code_native propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
# @show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
