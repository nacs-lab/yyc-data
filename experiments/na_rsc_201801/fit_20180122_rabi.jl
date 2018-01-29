#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Atomic: all_scatter_D
import NaCsCalc.Format: Unc, Sci
import NaCsCalc: Trap
using NaCsCalc.Utils: binomial_estimate, thread_rng, interactive
using NaCsData
using NaCsSim.DecayRabi: propagate_multistates, average_multistates, Γ_to_rates

using PyPlot
using NaCsPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180122_214615.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180122_111933.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)
data_b = NaCsData.select_count(params_b, logicals_b, selector)

const spec_a = OrderedDict(
    :xp1=>linspace(0, 135, 16),
    :yp1=>linspace(6, 90, 15),
    :zp1=>linspace(17, 255, 15),
    :x0=>linspace(3.5, 52.5, 15),
    :y0=>linspace(2.5, 37.5, 15),
    :z0=>linspace(7, 105, 15),
    :xf=>linspace(18, 20, 101),
    :yf=>linspace(18, 20, 101),
    :zf=>linspace(18.5, 19.1, 121)
)

const spec_b = OrderedDict(
    :xp1=>linspace(0, 270, 31),
    :yp1=>linspace(6, 180, 30),
    :zp1=>linspace(17, 510, 30),
    :x0=>linspace(3.5, 140, 40),
    :y0=>linspace(2.5, 100, 40),
    :z0=>linspace(7, 280, 40),
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)

const data_hot_xp1 = split_a[:xp1]
const data_hot_yp1 = [data_hot_xp1[1]; split_a[:yp1]]
const data_hot_zp1 = [data_hot_xp1[1]; split_a[:zp1]]
const data_hot_x0 = [data_hot_xp1[1]; split_a[:x0]]
const data_hot_y0 = [data_hot_xp1[1]; split_a[:y0]]
const data_hot_z0 = [data_hot_xp1[1]; split_a[:z0]]

const data_cold_xp1 = split_b[:xp1]
const data_cold_yp1 = [data_cold_xp1[1]; split_b[:yp1]]
const data_cold_zp1 = [data_cold_xp1[1]; split_b[:zp1]]
const data_cold_x0 = [data_cold_xp1[1]; split_b[:x0]]
const data_cold_y0 = [data_cold_xp1[1]; split_b[:y0]]
const data_cold_z0 = [data_cold_xp1[1]; split_b[:z0]]

const prefix = joinpath(@__DIR__, "imgs", "fit_20180122_rabi")

#### Fitting

# Scattering rates
# (detunings)
const δf1 = -75.0e9
const δf2 = -75.0e9 - 1.77e9
const rlof_f1 = (61.542e6 / (δf1 - 1.107266e9))^2
const rlof_f2 = (61.542e6 / (δf2 - 1.107266e9))^2
const rhif_f1 = (61.542e6 / (δf1 + 664.360e6))^2
const rhif_f2 = (61.542e6 / (δf2 + 664.360e6))^2

# polarizations of Raman beam, estimated
const rates_f1_coprop = all_scatter_D(true, 3, (0.5, 0.0, 0.5), rhif_f1, rlof_f1)
const rates_f1_up = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f1, rlof_f1)
const rates_f1_down = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f1, rlof_f1)
const rates_f2_coprop = all_scatter_D(true, 3, (0.25, 0.5, 0.25), rhif_f2, rlof_f2)
const rates_f2_counterop = all_scatter_D(true, 3, (0.1, 0.0, 0.9), rhif_f2, rlof_f2)

## TODO: check?
rates_f1_coprop .*= 4.46e8 / 3.18 * 5.38
rates_f2_coprop .*= 4.1e8 / 2.65 * 1.84
rates_f1_up .*= 1.05e9 / 4.24 * 6.31
rates_f1_down .*= 8.2e8 / 7.4 * 8.41
rates_f2_counterop .*= 3.25e9 / 0.25 * 2.69

const rates_rx = rates_f1_up + rates_f2_counterop
const rates_ry = rates_f1_down + rates_f2_counterop
const rates_az = rates_f1_up * (8 / 9) + rates_f2_coprop * 0.5
const rates_coprop = rates_f1_coprop + rates_f2_coprop

# Matrix elements
const m_Na = 23e-3 / 6.02e23
const η_rx = Trap.η(m_Na, 479e3, 2π / 589e-9) * √(2) * 0.96
const η_ry = Trap.η(m_Na, 492e3, 2π / 589e-9) * √(2) * 0.933
const η_az = Trap.η(m_Na, 85.7e3, 2π / 589e-9) * 0.67

const ns_az = 0:150
const meles_az_0 = Trap.sideband.(ns_az, ns_az, η_az)
const meles_az_p1 = Trap.sideband.(ns_az, ns_az .+ 1, η_az)

const ns_rx = 0:50
const meles_rx_0 = Trap.sideband.(ns_rx, ns_rx, η_rx)
const meles_rx_p1 = Trap.sideband.(ns_rx, ns_rx .+ 1, η_rx)

const ns_ry = 0:50
const meles_ry_0 = Trap.sideband.(ns_ry, ns_ry, η_ry)
const meles_ry_p1 = Trap.sideband.(ns_ry, ns_ry .+ 1, η_ry)

function f1_prob(Ωs, pΩ::AbstractArray, Γ::AbstractMatrix{T},
                 rates::AbstractVector{T}, tmax::T, atol=0.005, δΩ=T(0),
                 n::Integer=100000, rd=thread_rng(); offsetΩ=0) where T<:AbstractFloat
    nΩ = length(Ωs)
    nstates = length(rates)
    count = 0
    for i in 1:n
        r = rand(rd)
        j = 0
        @inbounds for _j in 1:nΩ
            j = _j
            r -= pΩ[j]
            if r < 0
                break
            end
        end
        Ω0 = Ωs[j]
        δ = δΩ * randn(rd) + offsetΩ
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
    (Ωs, pΩ, t, atol=0.005, δΩ=0; offsetΩ=0)->f1_prob(Ωs, pΩ, Γ32, rates32, t, atol, δΩ; offsetΩ=offsetΩ)
end

const f_rx = f1_prop_getter(rates_rx)
const f_ry = f1_prop_getter(rates_ry)
const f_az = f1_prop_getter(rates_az)

const img_survive = 0.95

function plot_f1(f::F, ts, _Ωs, _pΩ, δΩ=0, scale=img_survive; offset=0, offsetΩ=0, kws...) where F
    res = zeros(length(ts))
    Ωs = Float32.(_Ωs)
    pΩ = Float32.(_pΩ)
    @time Threads.@threads for i in 1:length(ts)
        res[i] = f(Ωs, pΩ, Float32(ts[i]), 0.002, Float32(δΩ); offsetΩ=offsetΩ) * scale + offset
    end
    plot(ts * 1e6, res; kws...)
end

function plot_f1_thermal(f, ts, Ωs, nbar, δΩ=0, scale=img_survive;
                         offset=0, offsetΩ=0, kws...)
    nstates = length(Ωs)
    ns = 0:(nstates - 1)
    pΩ = (nbar / (nbar + 1)).^ns ./ (nbar + 1)
    plot_f1(f, ts, Ωs, pΩ, δΩ, scale; offset=offset, offsetΩ=offsetΩ, kws...)
end

function plot_f1_thermal2(f, ts, Ωs, nbar, Ωs2, nbar2, δΩ=0, scale=img_survive;
                          offset=0, offsetΩ=0, kws...)
    nstates = length(Ωs)
    ns = 0:(nstates - 1)
    pΩ = (nbar / (nbar + 1)).^ns ./ (nbar + 1)
    nstates2 = length(Ωs2)
    ns2 = 0:(nstates2 - 1)
    pΩ2 = (nbar2 / (nbar2 + 1)).^ns2 ./ (nbar2 + 1)
    plot_f1(f, ts, Ωs * Ωs2', pΩ * pΩ2', δΩ, scale; offset=offset, offsetΩ=offsetΩ, kws...)
end

function diviation(f, data, Ωs, pΩ, δΩ=0, scale=1 / img_survive; offset=0, offsetΩ=0)
    params, ratios, uncs = NaCsData.get_values(data)
    perm = sortperm(params)
    params = params[perm] * 1e-6
    ratios = ratios[perm, 2] .* scale
    uncs = uncs[perm, 2] .* scale
    n = length(params)
    s = 0.0
    for i in 1:n
        s += ((f(Ωs, pΩ, Float32(params[i]), 0.001, δΩ;
                 offsetΩ=offsetΩ) - ratios[i] + offset) / uncs[i])^2
    end
    return s, n
end

#### Plotting

## X cold
# 97.5(20)

const τ_rx = 18.451e-6
const p_rx = [0.975, 0.10, 0.0]
const δΩ_rx = 15.54e3

# function div_rx_0(p0)
#     p = copy(p_rx)
#     p[1] = p0
#     s, n = diviation(f_rx, data_cold_x0, 2π / τ_rx * meles_rx_0[1:3], p, δΩ_rx)
#     return s / (n - 4)
# end

# function div_rx_p1(p0)
#     p = copy(p_rx)
#     p[1] = p0
#     s, n = diviation(f_rx, data_cold_xp1, 2π / τ_rx * meles_rx_p1[1:3], p, δΩ_rx)
#     return s / (n - 4)
# end

# @show div_rx_0(0.98)
# @show div_rx_p1(0.98)

# function diviation_rx(τ, p)
#     np = length(p)
#     d1, n1 = diviation(f_rx, data_cold_x0, 2π / τ * meles_rx_0[1:3], p, δΩ_rx)
#     d2, n2 = diviation(f_rx, data_cold_xp1, 2π / τ * meles_rx_p1[1:3], p, δΩ_rx)
#     return (d1 + d2) / (n1 + n2)
#     # return d1 / n1
# end
# function objective_rx(x)
#     r = diviation_rx(x[1] * 1e-6, [x[2:end]; 1.0])
#     @show x r
#     return r
# end
# init_params = [18.451, 0.975, 0.03]
# @show objective_rx(init_params)
# using Optim
# @show optimize(objective_rx, init_params)

# function diviation_hot_rx(τ, p, scale)
#     np = length(p)
#     d1, n1 = diviation(f_rx, data_hot_x0, 2π / τ * meles_rx_0, p, δΩ_rx, scale)
#     d2, n2 = diviation(f_rx, data_hot_xp1, 2π / τ * meles_rx_p1, p, δΩ_rx, scale)
#     return (d1 + d2) / (n1 + n2)
#     # return d1 / n1
# end
# function objective_hot_rx(x)
#     nstates = length(meles_rx_0)
#     ns = 0:(nstates - 1)
#     nbar = x[1]
#     scale = x[2]
#     pΩ = (nbar / (nbar + 1)).^ns ./ (nbar + 1)
#     r = diviation_hot_rx(x[3] * 1e-6, [pΩ; 1.0], scale)
#     @show x r
#     return r
# end
# init_params = [3.454, 1.043, 18.451]
# @show objective_hot_rx(init_params)
# using Optim
# @show optimize(objective_hot_rx, init_params)

# ps = linspace(0.9, 1.0, 41)
# plot(ps, div_rx_0.(ps), label="Carrier")
# plot(ps, div_rx_p1.(ps), label="Heating")
# legend()
# show()

const nbar_hot_rx = 3.459
const τ_hot_rx = 19.926e-6

figure()
ts_rx_p1 = linspace(0, 280e-6, 1001)
ts_rx_p1_hot = linspace(0, 150e-6, 501)
plot_f1(f_rx, ts_rx_p1, 2π / τ_rx * meles_rx_p1[1:3], p_rx, color="darkslateblue")
NaCsPlot.plot_survival_data(data_cold_xp1, fmt="C0o", label="Cold")
plot_f1_thermal(f_rx, ts_rx_p1_hot, 2π / τ_hot_rx * meles_rx_p1, nbar_hot_rx,
                δΩ_rx, 0.956528, color="C3")
NaCsPlot.plot_survival_data(data_hot_xp1, fmt="C1o", label="Hot")
grid()
ylim([0, 1])
title("X heating")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_xp1")

figure()
ts_rx_0 = linspace(0, 144e-6, 1001)
ts_rx_0_hot = linspace(0, 60e-6, 501)
plot_f1(f_rx, ts_rx_0, 2π / τ_rx * meles_rx_0[1:3], p_rx, color="darkslateblue")
NaCsPlot.plot_survival_data(data_cold_x0, fmt="C0o", label="Cold")
plot_f1_thermal(f_rx, ts_rx_0_hot, 2π / τ_hot_rx * meles_rx_0, nbar_hot_rx,
                δΩ_rx, 0.956528, color="C3")
NaCsPlot.plot_survival_data(data_hot_x0, fmt="C1o", label="Hot")
grid()
ylim([0, 1])
title("X carrier")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_x0")

## Y cold
# 95.0(30)

const τ_ry = 12.3e-6
const p_ry = [0.95, 0.10, 0.0]
const δΩ_ry = 15.54e3

# function div_ry_0(p0)
#     p = copy(p_ry)
#     p[1] = p0
#     s, n = diviation(f_ry, data_cold_y0, 2π / τ_ry * meles_ry_0[1:3], p, δΩ_ry)
#     return s / (n - 4)
# end

# function div_ry_p1(p0)
#     p = copy(p_ry)
#     p[1] = p0
#     s, n = diviation(f_ry, data_cold_yp1, 2π / τ_ry * meles_ry_p1[1:3], p, δΩ_ry)
#     return s / (n - 4)
# end

# @show div_ry_0(0.94)
# @show div_ry_p1(0.94)

# function diviation_ry(τ, p)
#     np = length(p)
#     d1, n1 = diviation(f_ry, data_cold_y0, 2π / τ * meles_ry_0[1:3], p, δΩ_ry)
#     d2, n2 = diviation(f_ry, data_cold_yp1, 2π / τ * meles_ry_p1[1:3], p, δΩ_ry)
#     return (d1 + d2) / (n1 + n2)
#     # return d1 / n1
# end
# function objective_ry(x)
#     r = diviation_ry(x[1] * 1e-6, [x[2:end]; 1.0])
#     @show x r
#     return r
# end
# init_params = [12.2, 0.975, 0.03] # [61.19, 14.49, 0.987, 0.01]
# @show objective_ry(init_params)
# using Optim
# @show optimize(objective_ry, init_params)

# function diviation_hot_ry(τ, p, scale)
#     np = length(p)
#     d1, n1 = diviation(f_ry, data_hot_y0, 2π / τ * meles_ry_0, p, δΩ_ry, scale)
#     d2, n2 = diviation(f_ry, data_hot_yp1, 2π / τ * meles_ry_p1, p, δΩ_ry, scale)
#     return (d1 + d2) / (n1 + n2)
#     # return d1 / n1
# end
# function objective_hot_ry(x)
#     nstates = length(meles_ry_0)
#     ns = 0:(nstates - 1)
#     nbar = x[1]
#     scale = x[2]
#     pΩ = (nbar / (nbar + 1)).^ns ./ (nbar + 1)
#     r = diviation_hot_ry(x[3] * 1e-6, [pΩ; 1.0], scale)
#     @show x r
#     return r
# end
# init_params = [3.454, 1.043, 15]
# @show objective_hot_ry(init_params)
# using Optim
# @show optimize(objective_hot_ry, init_params)

# ps = linspace(0.9, 1.0, 41)
# plot(ps, div_ry_0.(ps), label="Carrier")
# plot(ps, div_ry_p1.(ps), label="Heating")
# legend()
# show()

const nbar_hot_ry = 3.224
const τ_hot_ry = 12.347e-6

figure()
ts_ry_p1 = linspace(0, 180e-6, 1001)
ts_ry_p1_hot = linspace(0, 95e-6, 501)
plot_f1(f_ry, ts_ry_p1, 2π / τ_ry * meles_ry_p1[1:3], p_ry, color="darkslateblue")
NaCsPlot.plot_survival_data(data_cold_yp1, fmt="C0o", label="Cold")
plot_f1_thermal(f_ry, ts_ry_p1_hot, 2π / τ_hot_ry * meles_ry_p1, nbar_hot_ry,
                δΩ_ry, 0.97323, color="C3")
NaCsPlot.plot_survival_data(data_hot_yp1, fmt="C1o", label="Hot")
grid()
ylim([0, 1])
title("Y heating")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_yp1")

figure()
ts_ry_0 = linspace(0, 105e-6, 1001)
ts_ry_0_hot = linspace(0, 40e-6, 501)
plot_f1(f_ry, ts_ry_0, 2π / τ_ry * meles_ry_0[1:3], p_ry, color="darkslateblue")
NaCsPlot.plot_survival_data(data_cold_y0, fmt="C0o", label="Cold")
plot_f1_thermal(f_ry, ts_ry_0_hot, 2π / τ_hot_ry * meles_ry_0, nbar_hot_ry,
                δΩ_ry, 0.97323, color="C3")
NaCsPlot.plot_survival_data(data_hot_y0, fmt="C1o", label="Hot")
grid()
ylim([0, 1])
title("Y carrier")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_y0")

## Z cold
# 95.0(40)

const τ_az = 38.868e-6
const p_az = [0.95, 0.03, 0.02]
const δΩ_az = 16e3
const offsetΩ_az = 2e3

# function div_az_0(p0)
#     p = copy(p_az)
#     p[1] = p0
#     s, n = diviation(f_az, data_cold_z0, 2π / τ_az * (meles_az_0[1:3] * meles_ry_0[1:3]'),
#                      p * p_ry', δΩ_az)
#     return s / (n - 4)
# end

# function div_az_p1(p0)
#     p = copy(p_az)
#     p[1] = p0
#     s, n = diviation(f_az, data_cold_zp1, 2π / τ_az * (meles_az_p1[1:3] * meles_ry_0[1:3]'),
#                      p * p_ry', δΩ_az)
#     return s / (n - 4)
# end

# @show div_az_0(0.987)
# @show div_az_p1(0.987)

# function diviation_az(τ, p, δΩ, offsetΩ)
#     np = length(p)
#     d1, n1 = diviation(f_az, data_cold_z0, 2π / τ * meles_az_0[1:np] * meles_ry_0[1:3]',
#                        p * p_ry', δΩ, offsetΩ=offsetΩ)
#     d2, n2 = diviation(f_az, data_cold_zp1, 2π / τ * meles_az_p1[1:np] * meles_ry_0[1:3]',
#                        p * p_ry', δΩ, offsetΩ=offsetΩ)
#     return (d1 + d2) / (n1 + n2)
#     # return d1 / n1
# end
# function objective(x)
#     r = diviation_az(x[1] * 1e-6, [x[4:end]; 1.0], x[2] * 1e3, x[3] * 1e3)
#     @show x r
#     return r
# end
# init_params = [38.8, 15.4, 1.0, 0.97, 0.01]
# using Optim
# @show optimize(objective, init_params)

# function diviation_hot_az(τ, p, scale, δΩ, δΩ2)
#     np = length(p)
#     nstates_y = length(meles_ry_0)
#     ns_y = 0:(nstates_y - 1)
#     nbar_y = nbar_hot_ry
#     p_y = (nbar_y / (nbar_y + 1)).^ns_y ./ (nbar_y + 1)
#     d1, n1 = diviation(f_az, data_hot_y0, 2π / τ * meles_az_0 * meles_ry_0',
#                        p * p_y', δΩ, scale, offset=0.027)
#     d2, n2 = diviation(f_az, data_hot_yp1, 2π / τ * meles_az_p1 * meles_ry_0',
#                        p * p_y', δΩ2, scale, offset=0.027, offsetΩ=30e3)
#     return (d1 + d2) / (n1 + n2)
#     # return d1 / n1
# end
# function objective_hot_az(x)
#     nstates = length(meles_az_0)
#     ns = 0:(nstates - 1)
#     nbar = x[1]
#     scale = 1.043
#     δΩ = x[2]
#     τ = τ_az
#     δΩ2 = x[3]
#     pΩ = (nbar / (nbar + 1)).^ns ./ (nbar + 1)
#     r = diviation_hot_az(τ, [pΩ; 1.0], scale, δΩ, δΩ2)
#     @show x r
#     return r
# end
# init_params = [12.27, 30.77e3, 30.77e3 * 2.5]
# @show objective_hot_az(init_params)
# using Optim
# @show optimize(objective_hot_az, init_params)

# ps = linspace(0.9, 1.0, 41)
# plot(ps, div_az_0.(ps), label="Carrier")
# plot(ps, div_az_p1.(ps), label="Heating")
# legend()
# show()

const nbar_hot_az = 20.0
const τ_hot_az = τ_az
const δΩ_hot_az = 30.77e3 * 0.7
const δΩ_hot_az2 = 30.77e3 * 3

figure()
ts_az_p1 = linspace(0, 520e-6, 1001)
ts_az_p1_hot = linspace(0, 280e-6, 501)
plot_f1(f_az, ts_az_p1, 2π / τ_az * (meles_az_p1[1:3] * meles_ry_0[1:3]'),
        p_az * p_ry', δΩ_az, offsetΩ=offsetΩ_az, color="darkslateblue")
NaCsPlot.plot_survival_data(data_cold_zp1, fmt="C0o", label="Cold")
plot_f1_thermal2(f_az, ts_az_p1_hot, 2π / τ_hot_az * meles_az_p1, nbar_hot_az,
                 meles_ry_0, nbar_hot_ry, δΩ_hot_az, offset=0.027, color="C3")
NaCsPlot.plot_survival_data(data_hot_zp1, fmt="C1o", label="Hot")
grid()
ylim([0, 1])
title("Z heating")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_zp1")

figure()
ts_az_0 = linspace(0, 285e-6, 1001)
ts_az_0_hot = linspace(0, 120e-6, 501)
plot_f1(f_az, ts_az_0, 2π / τ_az * (meles_az_0[1:3] * meles_ry_0[1:3]'),
        p_az * p_ry', δΩ_az, offsetΩ=offsetΩ_az, color="darkslateblue")
NaCsPlot.plot_survival_data(data_cold_z0, fmt="C0o", label="Cold")
plot_f1_thermal2(f_az, ts_az_0_hot, 2π / τ_hot_az * meles_az_0, nbar_hot_az,
                 meles_ry_0, nbar_hot_ry, δΩ_hot_az2, offset=0.027, offsetΩ=30e3, color="C3")
NaCsPlot.plot_survival_data(data_hot_z0, fmt="C1o", label="Hot")
grid()
ylim([0, 1])
title("Z carrier")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_z0")

NaCsPlot.maybe_show()
