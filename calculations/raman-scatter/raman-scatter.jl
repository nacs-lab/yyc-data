#!/usr/bin/julia

# Compute Rabi flopping with the present of decay terms
# The Hamiltonian is assumed to be time independent and the Rabi drive is on-resonance

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))
using NaCsCalc.Utils: binomial_estimate
using NaCsCalc.Atomic: all_scatter_D
using NaCsSim.DecayRabi: propagate, average, average_multistates, Γ_to_rates
using PyPlot

const δf1 = -25.0e9
const δf2 = -25.0e9 - 1.77e9
const rlof_f1 = (61.542e6 / (δf1 - 1.107266e9))^2
const rlof_f2 = (61.542e6 / (δf2 - 1.107266e9))^2
const rhif_f1 = (61.542e6 / (δf1 + 664.360e6))^2
const rhif_f2 = (61.542e6 / (δf2 + 664.360e6))^2

const rates_f1_up = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f1, rlof_f1)
const rates_f1_down = rates_f1_up
const rates_f1_coprop = all_scatter_D(true, 3, (0.5, 0.0, 0.5), rhif_f1, rlof_f1)
const rates_f2_coprop = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f2, rlof_f2)
const rates_f2_counterop = all_scatter_D(true, 3, (0.1, 0.0, 0.9), rhif_f2, rlof_f2)

ts = linspace(0, 5e-5, 101)

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
f(ts, Γ, Ω, "blue")
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
# f2(Ω, Γ)

function f3(ts, Γ, Ω, color)
    res = Vector{Float64}(length(ts))
    unc = Vector{Float64}(length(ts))
    rates = Γ_to_rates(Γ)
    Ω32 = Float32(Ω)
    Γ32 = Float32.(Γ)
    rates32 = Float32.(rates)
    ts32 = Float32.(ts)
    res .= 0
    unc .= 0
    ntrial = 100000
    @time Threads.@threads for i in 1:length(ts)
        local n
        n = average_multistates(Ω32, 1, 2, Γ32, rates32, 1, ts32[i], ntrial)
        res[i], unc[i] = binomial_estimate(n[1] + n[3], ntrial)
    end
    errorbar(ts32, res, unc, fmt="^-", label="0", color=color)
end
Γ = [0 0e4 2e4
      1e4 0 0
      4e4 0 0] * 4
f3(ts, Γ, Ω, "cyan")
f3(ts, Γ, 0.0, "red")

ylim([0, 1])
legend()
grid()

show()
