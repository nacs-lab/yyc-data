#!/usr/bin/julia

# Compute Rabi flopping with the present of decay terms
# The Hamiltonian is assumed to be time independent and the Rabi drive is on-resonance

using NaCsSim.DecayRabi: propagate, average, Γ_to_rates
using PyPlot

δt = 1e-8
pts = 0:200:10000

function f(δt, pts, Γ, Ω, color)
    res = Vector{Float64}(length(pts))
    unc = Vector{Float64}(length(pts))
    rates = Γ_to_rates(Γ)
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
f(δt, pts, Γ, Ω, "blue")
Ω = 2π * 20e3
Γ = [0 1e4
      3e4 0] * 4
f(δt, pts, Γ, Ω, "red")
Ω = 2π * 2e3
Γ = [2e4 0
      0 3e4] * 4
f(δt, pts, Γ, Ω, "green")
Ω = 2π * 2e3
Γ = [0 1e4
      3e4 0] * 4
f(δt, pts, Γ, Ω, "orange")
# Γ = [0 3e4
#       3e4 0] * 2
# f(δt, pts, Γ, Ω, "red")
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
# f(δt, pts, Γ, Ω, "blue")
# ts = δt * pts
# function y(t, Ω)
#     Γ = 3e4 * 2
#     Ω′ = √(Ω^2 - Γ^2 / 4)
#     -exp(-Γ * t / 2) * (Γ / 2 / Ω′ * sin(Ω′ * t) + cos(Ω′ * t))
# end
# plot(ts, (1 .- y.(ts, Ω)) ./ 2, color="orange")
# Γ = [0 4e4
#       2e4 0] * 2
# f(δt, pts, Γ, Ω, "cyan")
# Γ = [0 6e4
#       0e4 0] * 2
# f(δt, pts, Γ, Ω, "orange")
# Γ = [0 3e4
#       3e4 0] * 2
# f(δt, pts, Γ, 0.1e3, "blue")

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
