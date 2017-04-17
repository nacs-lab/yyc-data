#!/usr/bin/julia

# Compute Rabi flopping with the present of decay terms
# The Hamiltonian is assumed to be time independent and the Rabi drive is on-resonance

using NaCsSim.DecayRabi: propagate, average, average_multistates, Γ_to_rates
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
        a, s = average(Ω32, Γ32, rates32, δt32 * pts[i], 100000)
        res[i] = a[1]
        unc[i] = s[1]
    end
    errorbar(pts * δt, res, unc, fmt="^-", label="0", color=color)
end
# Ω = 2π * 0.0001e3
# Γ = [0 1e4
#       3e4 0] * 4
# f(δt, pts, Γ, Ω, "red")
Ω = 2π * 0e3
Γ = [0 1e4
      3e4 0] * 4
f(δt, pts, Γ, Ω, "blue")
# Ω = 2π * 2e3
# Γ = [2e4 0
#       0 3e4] * 4
# f(δt, pts, Γ, Ω, "green")
# Ω = 2π * 2e3
# Γ = [0 1e4
#       3e4 0] * 4
# f(δt, pts, Γ, Ω, "orange")
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
# Ω = 2π * 400e3
# Γ = [3e4 0
#       0 3e4] * 2
# f(δt, pts, Γ, Ω, "blue")
# ts = δt * linspace(pts[1], pts[end], 1000)
# function y(t, Ω)
#     Γ = 3e4 * 2
#     Ω′ = √(Ω^2 - Γ^2 / 4)
#     -exp(-Γ * t / 2) * (Γ / 2 / Ω′ * sin(Ω′ * t) + cos(Ω′ * t))
# end
# plot(ts, (1 .- y.(ts, Ω)) ./ 2, color="orange")

function f3(δt, pts, Γ, Ω, color)
    res = Vector{Float64}(length(pts))
    unc = Vector{Float64}(length(pts))
    rates = Γ_to_rates(Γ)
    Ω32 = Float32(Ω)
    Γ32 = Float32.(Γ)
    rates32 = Float32.(rates)
    δt32 = Float32(δt)
    res .= 0
    unc .= 0
    @time Threads.@threads for i in 1:length(pts)
        local a, s
        a, s = average_multistates(Ω32, 1, 2, Γ32, rates32, 1, δt32 * pts[i], 100000)
        res[i] = a[1]
        unc[i] = s[1]
    end
    errorbar(pts * δt, res, unc, fmt="^-", label="0", color=color)
end
f3(δt, pts, Γ, Ω, "cyan")

ylim([0, 1])
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
# f2(Ω, Γ)

show()
