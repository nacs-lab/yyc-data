#!/usr/bin/julia

# Compute Rabi flopping with the present of decay terms
# The Hamiltonian is assumed to be time independent and the Rabi drive is on-resonance

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))
using NaCsCalc.Utils: binomial_estimate
using NaCsCalc.Atomic: all_scatter_D
import NaCsCalc: Trap
using NaCsSim.DecayRabi: propagate, average, average_multistates, Γ_to_rates
# using PyPlot

const δf1 = -25.0e9
const δf2 = -25.0e9 - 1.77e9
const rlof_f1 = (61.542e6 / (δf1 - 1.107266e9))^2
const rlof_f2 = (61.542e6 / (δf2 - 1.107266e9))^2
const rhif_f1 = (61.542e6 / (δf1 + 664.360e6))^2
const rhif_f2 = (61.542e6 / (δf2 + 664.360e6))^2

const rates_f1_coprop = all_scatter_D(true, 3, (0.5, 0.0, 0.5), rhif_f1, rlof_f1)
const rates_f1_up = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f1, rlof_f1)
const rates_f1_down = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f1, rlof_f1)
const rates_f2_coprop = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f2, rlof_f2)
const rates_f2_counterop = all_scatter_D(true, 3, (0.1, 0.0, 0.9), rhif_f2, rlof_f2)

rates_f1_coprop .*= 4.46e8
rates_f2_coprop .*= 4.1e8
rates_f1_up .*= 1.05e9
rates_f1_down .*= 8.2e8
rates_f2_counterop .*= 3.25e9

const rates_r3 = rates_f1_down + rates_f2_counterop
const rates_r2 = rates_f1_up + rates_f2_counterop
const rates_a1 = rates_f1_up * 0.52 / 5.3 + rates_f2_coprop * 0.356 / 0.77
const rates_coprop = rates_f1_coprop + rates_f2_coprop

const m_Na = 23e-3 / 6.02e23
const η_a1 = Trap.η(m_Na, 67e3, 2π / 589e-9) / √(2)
const η_r2 = Trap.η(m_Na, 420e3, 2π / 589e-9) * √(2)
const η_r3 = Trap.η(m_Na, 580e3, 2π / 589e-9) * √(2)

const ns_a1 = 0:30
const meles_a1_0 = Trap.sideband.(ns_a1, ns_a1, η_a1)
const meles_a1_p1 = Trap.sideband.(ns_a1, ns_a1 .+ 1, η_a1)

const ns_r2 = 0:10
const meles_r2_0 = Trap.sideband.(ns_r2, ns_r2, η_r2)
const meles_r2_p1 = Trap.sideband.(ns_r2, ns_r2 .+ 1, η_r2)

const ns_r3 = 0:10
const meles_r3_0 = Trap.sideband.(ns_r3, ns_r3, η_r3)
const meles_r3_p1 = Trap.sideband.(ns_r3, ns_r3 .+ 1, η_r3)

ts = linspace(0, 200e-6, 401)

function f(ts, Γ, Ω, color)
    res = Vector{Float64}(length(ts))
    unc = Vector{Float64}(length(ts))
    rates = Γ_to_rates(Γ)
    # res .= 0
    # unc .= 0
    # @time Threads.@threads for i in 1:length(ts)
    #     local a, s
    #     a, s = average(Ω, Γ, rates, ts[i], 10000)
    #     res[i] = a[1]
    #     unc[i] = s[1]
    # end
    # errorbar(ts, res, unc, fmt="-", label="0", color=color)
    Ω32 = Float32(Ω)
    Γ32 = Float32.(Γ)
    rates32 = Float32.(rates)
    ts32 = Float32.(ts)
    res .= 0
    unc .= 0
    @time Threads.@threads for i in 1:length(ts32)
        local a, s
        a, s = average(Ω32, Γ32, rates32, ts32[i], 100000)
        res[i] = a[1]
        unc[i] = s[1]
    end
    errorbar(ts32, res, unc, fmt="^-", label="0", color=color)
end
# Ω = 2π * 0.0001e3
# Γ = [0 1e4
#       3e4 0] * 4
# f(ts, Γ, Ω, "red")
Ω = 2π * 100e3
Γ = [0 1e4
      3e4 0] * 4
# f(ts, Γ, Ω, "blue")
# Ω = 2π * 2e3
# Γ = [2e4 0
#       0 3e4] * 4
# f(ts, Γ, Ω, "green")
# Ω = 2π * 2e3
# Γ = [0 1e4
#       3e4 0] * 4
# f(ts, Γ, Ω, "orange")
# Γ = [0 3e4
#       3e4 0] * 2
# f(ts, Γ, Ω, "red")
# plot_ts = linspace(ts[1], ts[end], 1000)
# function y(t, _Ω)
#     _Γ = 3e4 * 2
#     Ω = √(_Ω^2 + 2 * _Γ^2)
#     Γ = 3_Γ
#     Ω′ = √(Ω^2 - Γ^2 / 4)
#     -exp(-Γ * t / 2) * (Γ / 2 / Ω′ * sin(Ω′ * t) + cos(Ω′ * t))
# end
# plot(plot_ts, (1 .- y.(plot_ts, Ω)) ./ 2, color="orange")
# Ω = 2π * 400e3
# Γ = [3e4 0
#       0 3e4] * 2
# f(ts, Γ, Ω, "blue")
# function y(t, Ω)
#     Γ = 3e4 * 2
#     Ω′ = √(Ω^2 - Γ^2 / 4)
#     -exp(-Γ * t / 2) * (Γ / 2 / Ω′ * sin(Ω′ * t) + cos(Ω′ * t))
# end
# plot(plot_ts, (1 .- y.(plot_ts, Ω)) ./ 2, color="orange")

function f2(Ω, Γ)
    rates = Γ_to_rates(Γ)
    propagate(Ω, Γ, rates, 0.11e-3)
    @time Threads.@threads for i in 1:100000000
        propagate(Ω, Γ, rates, 0.11e-3)
    end
end
f2(Ω, Γ)

function f3(ts, Γ, Ωs, pΩ; kws...)
    res = Vector{Float64}(length(ts))
    unc = Vector{Float64}(length(ts))
    rates = Γ_to_rates(Γ)
    Ωs32 = Float32.(Ωs)
    Γ32 = Float32.(Γ)
    rates32 = Float32.(rates)
    ts32 = Float32.(ts)
    res .= 0
    unc .= 0
    ntrial = 100000
    @time Threads.@threads for i in 1:length(ts)
        local n
        n = average_multistates(Ωs32, pΩ, 1, 6, Γ32, rates32, 1, ts32[i], ntrial)
        res[i], unc[i] = binomial_estimate(n[6] + n[7] + n[8], ntrial)
    end
    errorbar(ts32 * 1e6, res, unc; kws...)
end

# const τ_r3 = 11.5e-6

# figure()
# ts_r3_0 = linspace(0, 80e-6, 201)
# f3(ts_r3_0, rates_r3, 2π / τ_r3 * meles_r3_0[1:3], [1.0, 0, 0],
#    fmt="-", color="blue", label="100%")
# f3(ts_r3_0, rates_r3, 2π / τ_r3 * meles_r3_0[1:3], [0.9, 0.1, 0],
#    fmt="-", color="red", label="90%")
# ylim([0, 1])
# xlim([ts_r3_0[1] * 1e6, ts_r3_0[end] * 1e6])
# title("Radial 3 carrier")
# legend()
# grid()

# figure()
# ts_r3_p1 = linspace(0, 180e-6, 201)
# f3(ts_r3_p1, rates_r3, 2π / τ_r3 * meles_r3_p1[1:3], [1.0, 0, 0],
#    fmt="-", color="cyan", label="100%")
# f3(ts_r3_p1, rates_r3, 2π / τ_r3 * meles_r3_p1[1:3], [0.9, 0.1, 0],
#    fmt="-", color="orange", label="90%")
# ylim([0, 1])
# xlim([ts_r3_p1[1] * 1e6, ts_r3_p1[end] * 1e6])
# title("Radial 3 heating")
# legend()
# grid()

# const τ_r2 = 11.55e-6

# figure()
# ts_r2_0 = linspace(0, 80e-6, 201)
# f3(ts_r2_0, rates_r2, 2π / τ_r2 * meles_r2_0[1:3], [1.0, 0, 0],
#    fmt="-", color="blue", label="100%")
# f3(ts_r2_0, rates_r2, 2π / τ_r2 * meles_r2_0[1:3], [0.9, 0.05, 0.05],
#    fmt="-", color="red", label="90%")
# ylim([0, 1])
# xlim([ts_r2_0[1] * 1e6, ts_r2_0[end] * 1e6])
# title("Radial 2 carrier")
# legend()
# grid()

# figure()
# ts_r2_p1 = linspace(0, 180e-6, 201)
# f3(ts_r2_p1, rates_r2, 2π / τ_r2 * meles_r2_p1[1:3], [1.0, 0, 0],
#    fmt="-", color="cyan", label="100%")
# f3(ts_r2_p1, rates_r2, 2π / τ_r2 * meles_r2_p1[1:3], [0.9, 0.05, 0.05],
#    fmt="-", color="orange", label="90%")
# ylim([0, 1])
# xlim([ts_r2_p1[1] * 1e6, ts_r2_p1[end] * 1e6])
# title("Radial 2 heating")
# legend()
# grid()

# const τ_a1 = 61.1e-6

# figure()
# ts_a1_0 = linspace(0, 300e-6, 201)
# f3(ts_a1_0, rates_a1, 2π / τ_a1 * (meles_a1_0[1:3] * meles_r3_0[1:2]'),
#    [1.0, 0.0, 0.0] * [0.9, 0.1]',
#    fmt="-", color="blue", label="100%")
# f3(ts_a1_0, rates_a1, 2π / τ_a1 * (meles_a1_0[1:3] * meles_r3_0[1:2]'),
#    [0.9, 0.1, 0.0] * [0.9, 0.1]',
#    fmt="-", color="red", label="90%")
# ylim([0, 1])
# xlim([ts_a1_0[1] * 1e6, ts_a1_0[end] * 1e6])
# title("Axial carrier")
# legend()
# grid()

# figure()
# ts_a1_p1 = linspace(0, 450e-6, 201)
# f3(ts_a1_p1, rates_a1, 2π / τ_a1 * (meles_a1_p1[1:3] * meles_r3_0[1:2]'),
#    [1.0, 0.0, 0.0] * [0.9, 0.1]',
#    fmt="-", color="blue", label="100%")
# f3(ts_a1_p1, rates_a1, 2π / τ_a1 * (meles_a1_p1[1:3] * meles_r3_0[1:2]'),
#    [0.9, 0.1, 0.0] * [0.9, 0.1]',
#    fmt="-", color="orange", label="90%")
# ylim([0, 1])
# xlim([ts_a1_p1[1] * 1e6, ts_a1_p1[end] * 1e6])
# title("Axial heating")
# legend()
# grid()

# show()
