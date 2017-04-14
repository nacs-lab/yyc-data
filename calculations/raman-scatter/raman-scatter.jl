#!/usr/bin/julia

# Compute Rabi flopping with the present of decay terms
# The Hamiltonian is assumed to be time independent and the Rabi drive is on-resonance

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
    ΓΩ′::T # Or ΓΔ′
    Ω²pΓ²mΔ²::T
    p2max::T
    overdamp::Bool
end

function RabiDecayParams(Ω::T, Γ₁::T, Γ₂::T) where T
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
    Ω²pΓ²mΔ² = Ω² + Γ² - Δ²
    p2max = Ω^2 * Γ₂ / 2 / Γ / Ω²pΓ²mΔ²
    return RabiDecayParams{T}(Γ₁, Γ₂, Δ, Ω, Γ, Δ², Ω², Γ², Ω′², Ω′, Δ * Ω′, Γ * Ω′, Ω²pΓ²mΔ²,
                              p2max, overdamp)
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
function propagate_2states_underdamp{T}(params::RabiDecayParams{T}, tmax, rd)
    r = T(rand(rd)) * params.Ω′²
    # Now find the t for which `ψ^2(t) * Ω′² = r`.
    # First check if `ψ^2(tmax) * Ω′² > r`
    t::T = tmax
    Ω′t = params.Ω′ * t
    @fastmath sinΩ′t = sin(Ω′t)
    @fastmath cosΩ′t = cos(Ω′t)
    @fastmath expΓt = exp(-params.Γ * t)
    ψ² = -expΓt * muladd(params.ΔΩ′, sinΩ′t, muladd(params.Δ², cosΩ′t, -params.Ω²))
    if ψ² > r
        # No decay happened, return the wave function at tmax
        Ω′t_2 = Ω′t / 2
        @fastmath sinΩ′t_2 = sin(Ω′t_2)
        @fastmath cosΩ′t_2 = cos(Ω′t_2)
        ψ1 = muladd(params.Ω′, cosΩ′t_2, -params.Δ * sinΩ′t_2)
        ψ2 = params.Ω * sinΩ′t_2
        @fastmath factor = sqrt(expΓt / ψ²)
        return (t, 0, (ψ1 * factor, ψ2 * factor))
    end
    # Tolerance
    yδ = max(T(2e-7), eps(T) * 10) * params.Ω′²
    # Find the root using a combination of newton's method and bisecting.
    # The bisection is to avoid newton's method overshotting since the function we want
    # to solve is known to have partial oscillation.
    tlo::T = 0
    thi::T = t
    if r - ψ² <= yδ
        @goto ret
    end
    if params.Ω′² - r <= yδ
        return (tlo, 1, (one(T), zero(T)))
    end
    if params.Ω′ > 0
        tthresh = min(T(2 / params.Ω′), T(tmax))
    else
        tthresh = T(tmax)
    end
    c1 = params.Γ * params.Ω²
    c2 = params.ΔΩ′ * params.Γ₁
    c3 = params.Δ * muladd(params.Γ, -params.Δ, params.Ω′²)
    while thi - tlo > 5 * eps(T) * thi
        if thi - tlo < tthresh
            # Try Newton's method
            diff = -expΓt * muladd(c3, cosΩ′t, -muladd(c2, sinΩ′t, - c1))
            δ1 = (ψ² - r)
            t2 = t - δ1 / diff
            if tlo < t2 < thi
                δt = thi - tlo
                t = t2
                Ω′t = params.Ω′ * t
                @fastmath sinΩ′t = sin(Ω′t)
                @fastmath cosΩ′t = cos(Ω′t)
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
        t = (thi + tlo) / 2
        Ω′t = params.Ω′ * t
        @fastmath sinΩ′t = sin(Ω′t)
        @fastmath cosΩ′t = cos(Ω′t)
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
    p2 = params.p2max * (1 - expΓt / params.Ω′² * muladd(params.ΓΩ′, sinΩ′t,
                                                            -(muladd(params.Γ², cosΩ′t,
                                                                     -params.Ω²pΓ²mΔ²))))
    ptotal = 1 - ψ² / params.Ω′²
    return (t, ptotal * rand(rd) < p2 ? 2 : 1, (one(T), zero(T)))
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
    @fastmath sinΩt_2 = sin(Ω * δt / 2)
    @fastmath cosΩt_2 = cos(Ω * δt / 2)
    # Full step
    @fastmath sinΩt = sin(Ω * δt)
    @fastmath cosΩt = cos(Ω * δt)
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
            v = ϕ[j] * (1 - rates[j] * δt) * factor
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


# ϕ = [1.0, 0.0]
# δt = 1e-7
# Γ = [0 3e4
#       2e4 0]
# i1 = 1
# i2 = 2
# Ω = 2π * 400e3

# using PyPlot
# pts = 0:10:1000
# res = Vector{Float64}(length(pts))
# unc = Vector{Float64}(length(pts))

# function f(pts, ϕ, Γ, i1, i2, Ω, res, unc)
#     len = length(ϕ)
#     rates = Vector{eltype(Γ)}(len)
#     for i in 1:len
#         local s = zero(eltype(Γ))
#         for j in 1:len
#             s += Γ[j, i]
#         end
#         rates[i] = s
#     end
#     rds = [MersenneTwister(0) for i in 1:Threads.nthreads()]
#     # @show average(ϕ, 1.1e-7, 0, Γ, i1, i2, Ω, 10000, rates, rds[Threads.threadid()])
#     @time Threads.@threads for i in 1:length(pts)
#         local a, s
#         a, s = average(ϕ, 1.1e-7, pts[i], Γ, i1, i2, Ω, 10000, rates, rds[Threads.threadid()])
#         res[i] = a[1]
#         unc[i] = s[1]
#     end
#     errorbar(pts, res, unc)
#     @time Threads.@threads for i in 1:length(pts)
#         local a, s
#         a, s = average(ϕ, 1.1e-7 * 5, pts[i] ÷ 5, Γ, i1, i2, Ω, 10000, rates, rds[Threads.threadid()])
#         res[i] = a[1]
#         unc[i] = s[1]
#     end
#     errorbar(pts, res, unc)
# end
# f(pts, ϕ, Γ, i1, i2, Ω, res, unc)
# # # for i in 1:npts
# # #     a, s = average(ϕ, 1.1e-7, 10 * (i - 1), Γ, i1, i2, 0, 10000)
# # #     res[i] = a[1]
# # #     unc[i] = s[1]
# # # end
# # # errorbar(1:npts, res, unc)

# grid()
# show()

# # @time @show average(ϕ, 5.5e-7, 200, Γ, i1, i2, 0, 1)
# # @time @show average(ϕ, 2.75e-7, 40_000, Γ, i1, i2, Ω, 100)
# # @time @show average(ϕ, 1.1e-7, 100_000, Γ, i1, i2, Ω, 100)

params = RabiDecayParams(2.0, 0.5, 0.4)

@code_native propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
@show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
@show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
@show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
@show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
@show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
@show propagate_2states_underdamp(params, 100, Base.Random.GLOBAL_RNG)
