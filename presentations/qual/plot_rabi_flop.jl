#!/usr/bin/julia

# Compute Rabi flopping with the present of decay terms
# The Hamiltonian is assumed to be time independent and the Rabi drive is on-resonance

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))
using NaCsCalc.Utils: binomial_estimate, thread_rng, interactive
using NaCsCalc.Atomic: all_scatter_D
import NaCsCalc: Trap
using NaCsSim.DecayRabi: propagate_multistates, average_multistates, Γ_to_rates
using NaCsData
using DataStructures
using PyPlot
matplotlib["rcParams"][:update](Dict("font.size" => 20))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

# Scattering rates
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

# Matrix elements
const m_Na = 23e-3 / 6.02e23
const η_a1 = Trap.η(m_Na, 68.8e3, 2π / 589e-9) * 0.67
const η_r2 = Trap.η(m_Na, 430e3, 2π / 589e-9) * √(2)
const η_r3 = Trap.η(m_Na, 589.5e3, 2π / 589e-9) * √(2)

const ns_a1 = 0:150
const meles_a1_0 = Trap.sideband.(ns_a1, ns_a1, η_a1)
const meles_a1_p1 = Trap.sideband.(ns_a1, ns_a1 .+ 1, η_a1)

const ns_r2 = 0:50
const meles_r2_0 = Trap.sideband.(ns_r2, ns_r2, η_r2)
const meles_r2_p1 = Trap.sideband.(ns_r2, ns_r2 .+ 1, η_r2)

const ns_r3 = 0:50
const meles_r3_0 = Trap.sideband.(ns_r3, ns_r3, η_r3)
const meles_r3_p1 = Trap.sideband.(ns_r3, ns_r3 .+ 1, η_r3)

const iname_a = joinpath(@__DIR__, "../../experiments/misc/data/data_20170404_133229.csv")
const iname_b = joinpath(@__DIR__, "../../experiments/misc/data/data_20170405_183620.csv")
const iname_c = joinpath(@__DIR__, "../../experiments/misc/data/data_20170409_141913.csv")
const iname_d = joinpath(@__DIR__, "../../experiments/misc/data/data_20170409_162137.csv")
const iname_e = joinpath(@__DIR__, "../../experiments/misc/data/data_20170409_184250.csv")
const iname_f = joinpath(@__DIR__, "../../experiments/misc/data/data_20170409_233639.csv")

# Before Carrier
data_a = NaCsData.load_count_csv(iname_a)
# Before +1
data_b = NaCsData.load_count_csv(iname_b)
# After a1 +1
data_c = NaCsData.load_count_csv(iname_c)
# After a1 0
data_d = NaCsData.load_count_csv(iname_d)
# After r2 0, +1 (<0 is carrier)
data_e = NaCsData.load_count_csv(iname_e)
# After r3 0, +1
data_f = NaCsData.load_count_csv(iname_f)

# Before Carrier
spec_a = OrderedDict(
    :unused=>(linspace(-17.53, -17.745, 11),
              linspace(-17.14, -17.44, 11),
              linspace(-18.40, -17.90, 51),
              linspace(-18.40, -18.35, 11)),
    :after=>(0:2:100, 2:2:100, 10:10:400),
    :before=>(2:2:30, 2:2:30, 10:10:140)
)
# Before +1
spec_b = OrderedDict(
    :after=>(0:7:280, 8:8:240, 25:25:450),
    :before=>(7:7:140, 8:8:160, 25:25:400)
)
# After r3 0, +1
spec_f = (0:6:180, 2:2:80)

data_a = NaCsData.split_data(data_a, spec_a)
data_b = NaCsData.split_data(data_b, spec_b)
data_f = NaCsData.split_data(data_f, spec_f)

function plot_data(data, scale=1; kws...)
    params, ratios, uncs = NaCsData.get_values(data)
    perm = sortperm(params)
    params = params[perm]
    ratios = ratios[perm, 2] .* scale
    uncs = uncs[perm, 2] .* scale
    errorbar(params, ratios, uncs; kws...)
end

function maybe_save(name)
    if !interactive()
        dir = dirname(name)
        if !isempty(dir)
            mkpath(dir, 0o755)
        end
        savefig("$name.pdf"; bbox_inches="tight", transparent=true)
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg", bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if interactive()
        show()
    end
end

const prefix = joinpath(@__DIR__, "imgs/rabi_flop")

data_before_r2_0 = [data_a[:after][1][1]; data_a[:before][1]]
data_before_r3_0 = [data_a[:after][1][1]; data_a[:before][2]]
data_before_a1_0 = [data_a[:after][1][1]; data_a[:before][3]]

data_before_r2_p1 = [data_b[:after][1][1]; data_b[:before][1]]
data_before_r3_p1 = [data_b[:after][1][1]; data_b[:before][2]]
data_before_a1_p1 = [data_b[:after][1][1]; data_b[:before][3]]

data_after_r2_0 = NaCsData.map_params((i, v)->-v * 1e6, data_e[data_e.params .<= 0])
data_after_r2_p1 = NaCsData.map_params((i, v)->v * 1e6, data_e[data_e.params .>= 0])

data_after_r3_0 = [data_f[1][1]; data_f[2]]
data_after_r3_p1 = data_f[1]

data_after_a1_0 = NaCsData.map_params((i, v)->v * 1e6, data_d)
data_after_a1_p1 = NaCsData.map_params((i, v)->v * 1e6, data_c)

function f1_prob(Ωs, pΩ::AbstractArray, Γ::AbstractMatrix{T},
                 rates::AbstractVector{T}, tmax::T, atol=0.005, δΩ=T(0),
                 n::Integer=100000, rd=thread_rng()) where T<:AbstractFloat
    nΩ = length(Ωs)
    nstates = length(rates)
    count = 0
    for i in 1:n
        r = rand(rd)
        j = 0
        @inbounds for j in 1:nΩ
            r -= pΩ[j]
            if r < 0
                break
            end
        end
        Ω0 = Ωs[j]
        δ = δΩ * randn(rd)
        Ω = T(sqrt(Ω0^2 + δ^2))
        i_final = propagate_multistates(Ω, 1, 6, Γ, rates, 1, tmax, rd)
        if i_final > 5
            if rand(rd) < Ω0^2 / Ω^2
                # count F1
                count += 1
            end
        end
        if i % 256 == 0
            r, s = binomial_estimate(count, i)
            if s < atol
                return r
            end
        end
    end
    return binomial_estimate(count, n)[1]
end

function f1_prop_getter(Γ)
    Γ32 = Float32.(Γ)
    rates32 = Γ_to_rates(Γ32)
    (Ωs, pΩ, t, atol=0.005, δΩ=0)->f1_prob(Ωs, pΩ, Γ32, rates32, t, atol, δΩ)
end

const f_r3 = f1_prop_getter(rates_r3)
const f_r2 = f1_prop_getter(rates_r2)
const f_a1 = f1_prop_getter(rates_a1)

function plot_f1(f, ts, Ωs, pΩ, δΩ=0, scale=0.85; kws...)
    res = zeros(length(ts))
    @time Threads.@threads for i in 1:length(ts)
        res[i] = f(Ωs, pΩ, Float32(ts[i]), 0.002, δΩ) * scale
    end
    plot(ts * 1e6, res; kws...)
end

function plot_f1_thermal(f, ts, Ωs, nbar, δΩ=0, scale=0.92; kws...)
    nstates = length(Ωs)
    ns = 0:(nstates - 1)
    pΩ = (nbar / (nbar + 1)).^ns ./ (nbar + 1)
    plot_f1(f, ts, Ωs, pΩ, δΩ, scale; kws...)
end

function plot_f1_thermal2(f, ts, Ωs, nbar, Ωs2, nbar2, δΩ=0, scale=0.92; kws...)
    nstates = length(Ωs)
    ns = 0:(nstates - 1)
    pΩ = (nbar / (nbar + 1)).^ns ./ (nbar + 1)
    nstates2 = length(Ωs2)
    ns2 = 0:(nstates2 - 1)
    pΩ2 = (nbar2 / (nbar2 + 1)).^ns2 ./ (nbar2 + 1)
    plot_f1(f, ts, Ωs * Ωs2', pΩ * pΩ2', δΩ, scale; kws...)
end

# function diviation(f, data, Ωs, pΩ, scale=1 / 0.85, δΩ=0)
#     params, ratios, uncs = NaCsData.get_values(data)
#     perm = sortperm(params)
#     params = params[perm] * 1e-6
#     ratios = ratios[perm, 2] .* scale
#     uncs = uncs[perm, 2] .* scale
#     n = length(params)
#     s = 0.0
#     for i in 1:n
#         s += ((f(Ωs, pΩ, Float32(params[i]), 0.001, δΩ) - ratios[i]) / uncs[i])^2
#     end
#     return s, n
# end

# r3 0.90±0.05

const τ_r3 = 11.445e-6
const p_r3 = [0.93, 0.057, 0.013]
const δΩ_r3 = 0

const τ_r3_init = 12.8e-6
const nbar_init_r3 = 2.5
const δΩ_r3_init = 15e3

figure()
ts_r3_0 = linspace(0, 80e-6, 256)
plot_data(data_after_r3_0, 1, fmt="C0s")
plot_f1(f_r3, ts_r3_0, 2π / τ_r3 * meles_r3_0[1:3], p_r3, color="C0")
ts_r3_0_init = linspace(0, 30e-6, 256)
plot_data(data_before_r3_0, 1, fmt="C3o")
plot_f1_thermal(f_r3, ts_r3_0_init, 2π / τ_r3_init * meles_r3_0, nbar_init_r3, δΩ_r3_init,
                color="C3")
axhline(0.85, c="C0", ls="-.")
ylim([0, 1])
xlim([ts_r3_0[1] * 1e6, ts_r3_0[end] * 1e6])
xlabel("Time (\$\\mu s\$)")
ylabel("F=1 population")
grid()
maybe_save("$(prefix)_r3_0")

figure()
ts_r3_p1 = linspace(0, 180e-6, 256)
plot_data(data_after_r3_p1, 1, fmt="C0s")
plot_f1(f_r3, ts_r3_p1, 2π / τ_r3 * meles_r3_p1[1:3], p_r3, color="C0")
ts_r3_p1_init = linspace(0, 160e-6, 256)
plot_data(data_before_r3_p1, 1, fmt="C3o")
plot_f1_thermal(f_r3, ts_r3_p1_init, 2π / τ_r3_init * meles_r3_p1, nbar_init_r3, δΩ_r3_init, color="C3")
axhline(0.85, c="C0", ls="-.")
ylim([0, 1])
xlim([ts_r3_p1[1] * 1e6, ts_r3_p1[end] * 1e6])
xlabel("Time (\$\\mu s\$)")
ylabel("F=1 population")
grid()
maybe_save("$(prefix)_r3_p1")

# r2 0.90±0.06

const τ_r2 = 11.608e-6
const p_r2 = [0.896, 0.048, 0.056]
const δΩ_r2 = 0

const τ_r2_init = 13.5e-6
const nbar_init_r2 = 3.6
const δΩ_r2_init = 15e3

figure()
ts_r2_0 = linspace(0, 80e-6, 256)
plot_data(data_after_r2_0, 1, fmt="C0s")
plot_f1(f_r2, ts_r2_0, 2π / τ_r2 * meles_r2_0[1:3], p_r2, color="C0")
ts_r2_0_init = linspace(0, 30e-6, 256)
plot_data(data_before_r2_0, 1, fmt="C3o")
plot_f1_thermal(f_r2, ts_r2_0_init, 2π / τ_r2_init * meles_r2_0, nbar_init_r2,
                δΩ_r2_init, color="C3")
axhline(0.85, c="C0", ls="-.")
ylim([0, 1])
xlim([ts_r2_0[1] * 1e6, ts_r2_0[end] * 1e6])
xlabel("Time (\$\\mu s\$)")
ylabel("F=1 population")
grid()
maybe_save("$(prefix)_r2_0")

figure()
ts_r2_p1 = linspace(0, 180e-6, 256)
plot_data(data_after_r2_p1, 1, fmt="C0s")
plot_f1(f_r2, ts_r2_p1, 2π / τ_r2 * meles_r2_p1[1:3], p_r2, color="C0")
ts_r2_p1_init = linspace(0, 140e-6, 256)
plot_data(data_before_r2_p1, 1, fmt="C3o")
plot_f1_thermal(f_r2, ts_r2_p1_init, 2π / τ_r2_init * meles_r2_p1, nbar_init_r2,
                δΩ_r2_init, color="C3")
axhline(0.85, c="C0", ls="-.")
ylim([0, 1])
xlim([ts_r2_p1[1] * 1e6, ts_r2_p1[end] * 1e6])
xlabel("Time (\$\\mu s\$)")
ylabel("F=1 population")
grid()
maybe_save("$(prefix)_r2_p1")

# a1 0.91

const τ_a1 = 60.3e-6
const p_a1 = [0.92, 0.05, 0.02]
const δΩ_a1 = 15e3

const τ_a1_init = 80e-6
const nbar_init_a1 = 10.0
const δΩ_a1_init_0 = 25e3
const δΩ_a1_init_p1 = 50e3

figure()
ts_a1_0 = linspace(0, 300e-6, 256)
plot_data(data_after_a1_0, 1, fmt="C0s")
plot_f1(f_a1, ts_a1_0, 2π / τ_a1 * (meles_a1_0[1:3] * meles_r3_0[1:3]'), p_a1 * p_r3', color="C0")
ts_a1_0_init = linspace(0, 150e-6, 256)
plot_data(data_before_a1_0, 1, fmt="C3o")
plot_f1_thermal2(f_a1, ts_a1_0_init, 2π / τ_a1_init * meles_a1_0, nbar_init_a1,
                 meles_r3_0, nbar_init_r3, δΩ_a1_init_0, color="C3")
axhline(0.85, c="C0", ls="-.")
ylim([0, 1])
xlim([ts_a1_0[1] * 1e6, ts_a1_0[end] * 1e6])
xlabel("Time (\$\\mu s\$)")
ylabel("F=1 population")
grid()
maybe_save("$(prefix)_a1_0")

figure()
ts_a1_p1 = linspace(0, 450e-6, 256)
plot_data(data_after_a1_p1, 1, fmt="C0s")
plot_f1(f_a1, ts_a1_p1, 2π / τ_a1 * (meles_a1_p1[1:3] * meles_r3_0[1:3]'),
        p_a1 * p_r3', 0, color="C0")
plot_f1(f_a1, ts_a1_p1, 2π / τ_a1 * (meles_a1_p1[1:3] * meles_r3_0[1:3]'),
        p_a1 * p_r3', δΩ_a1, color="C2", ls="-.")
ts_a1_p1_init = linspace(0, 400e-6, 256)
plot_data(data_before_a1_p1, 1, fmt="C3o")
plot_f1_thermal2(f_a1, ts_a1_p1_init, 2π / τ_a1_init * meles_a1_p1, nbar_init_a1,
                 meles_r3_0, nbar_init_r3, δΩ_a1_init_p1, color="C3")
axhline(0.85, c="C0", ls="-.")
ylim([0, 1])
xlim([ts_a1_p1[1] * 1e6, ts_a1_p1[end] * 1e6])
xlabel("Time (\$\\mu s\$)")
ylabel("F=1 population")
grid()
maybe_save("$(prefix)_a1_p1")

maybe_show()
