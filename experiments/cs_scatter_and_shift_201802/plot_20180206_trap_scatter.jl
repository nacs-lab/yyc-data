#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Atomic
import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180206_150235.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180206_164731.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[2, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[4, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)
data_b = NaCsData.select_count(params_b, logicals_b, selector)

const spec_a = OrderedDict(
    :f3=>[0.001, 50, 100, 200, 400, 800],
    :f44=>[0.001, 50, 100, 200, 400, 800],
    :all=>[0.001, 50, 100, 200, 400, 800],
)

const spec_b = OrderedDict(
    :f3=>[1600.0, 2400.0],
    :f44=>[1600.0, 2400.0],
    :all=>[1600.0, 2400.0],
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)

const data_f3 = [split_a[:f3]; split_b[:f3]]
const data_f44 = [split_a[:f44]; split_b[:f44]]
const data_all = [split_a[:all]; split_b[:all]]

const prefix = joinpath(@__DIR__, "imgs", "data_20180206_trap_scatter")

function get_survival_data(data, z=1.0)
    params, ratios, uncs = NaCsData.get_values(data, z)
    perm = sortperm(params)
    params = params[perm]
    return params, ratios[perm, 2], uncs[perm, 2]
end

function get_survival_ratios(base, sub, z=1.0)
    params_b, ratios_b, uncs_b = get_survival_data(base, z)
    params_s, ratios_s, uncs_s = get_survival_data(sub, z)
    params_b, ratios_s ./ ratios_b, sqrt.((uncs_s ./ ratios_b).^2 .+
                                          (uncs_b ./ ratios_b.^2 .* ratios_s).^2)
end

function plot_survival_ratios(base, sub, z=1.0; kws...)
    params, ratios, uncs = get_survival_ratios(base, sub, z)
    errorbar(params, ratios, uncs; kws...)
end

const cs_atom = Atomic.Alkali(7, 335.116048807e12, 351.72571850e12, 28.743e6, 32.889e6,
                              (4.021776399375e9, -5.170855370625e9),
                              (510.860e6, -656.820e6),
                              (263.8906e6, 12.79851e6, -188.4885e6, -339.7128e6))

function all_scatters(atom, Ω, freq, pol)
    # Order: F high->low; mF low->high
    idx_to_state = function (idx)
        local Fx2, mFx2, r
        if idx > atom.Ix2 + 2
            # low F
            idx -= atom.Ix2 + 2
            Fx2 = atom.Ix2 - 1
            mFx2 = idx * 2 - atom.Ix2 - 1
        else
            # high F
            Fx2 = atom.Ix2 + 1
            mFx2 = idx * 2 - atom.Ix2 - 3
        end
        return Fx2, mFx2
    end
    nstates = (atom.Ix2 + 1) * 2
    rates = Matrix{Float64}(nstates, nstates)
    @inbounds for i in 1:nstates
        F1x2, mF1x2 = idx_to_state(i)
        for j in 1:nstates
            F2x2, mF2x2 = idx_to_state(j)
            rates[j, i] = Atomic.get_scatter(atom, Ω, freq, F1x2, mF1x2, F2x2, mF2x2, pol)
        end
    end
    return rates
end

const rates_trap = all_scatters(cs_atom, 1e10, 299792458 / 976e-9, (0, 1, 0))

function rates_to_A(rates)
    nx, ny = size(rates)
    A = Matrix{Float64}(nx, ny)
    @inbounds for i in 1:nx
        s = 0.0
        for j in 1:ny
            r = rates[j, i]
            A[j, i] = r
            s += r
        end
        A[i, i] -= s
    end
    return A
end

function propagate(A, init, t, idxs)
    res = exp(A * t) * init # expm
    sum(getindex.((res,), idxs))
end

function gen_propagate(rates, init, idxs)
    A = rates_to_A(rates)
    t->propagate(A, init, t, idxs)
end

function gen_model(rates, init, idxs)
    propagate_f = gen_propagate(rates, init, idxs)
    (x, p) -> p[1] .* propagate_f.(x .* p[2])
end

function gen_model2(rates, init, idxs)
    propagate_f = gen_propagate(rates, init, idxs)
    (x, p) -> propagate_f.(x .* p[1])
end

function loss_model(x, p)
    return p[1] .* exp.(-x ./ p[2])
end

const init_mf = [0, 0, 0, 0, 0, 0, 0, 0, 1.0,
                 0, 0, 0, 0, 0, 0, 0]
const model_f3 = gen_model2(rates_trap, init_mf, (10, 11, 12, 13, 14, 15, 16))
# const model_f44 = gen_model(rates_trap, init_mf, (9,))

function fit_survival_ratios(model, base, sub, p0; plotx=nothing, use_unc=false, plot_scale=1.1)
    if use_unc
        params, ratios, uncs = get_survival_ratios(base, sub)
    else
        params, ratios, uncs = get_survival_ratios(base, sub, 0.0)
    end
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        plotx = linspace(mid - span * plot_scale / 2, mid + span * plot_scale / 2, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios, 1 ./ uncs.^2, p0)
    else
        fit = curve_fit(model, params, ratios, p0)
    end
    return (param=fit.param, unc=estimate_errors(fit),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

function fit_survival(model, data, p0; plotx=nothing, use_unc=false, plot_scale=1.1)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
    end
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        plotx = linspace(mid - span * plot_scale / 2, mid + span * plot_scale / 2, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios[:, 2], 1 ./ uncs[:, 2].^2, p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
    return (param=fit.param, unc=estimate_errors(fit),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

const fit_f3 = fit_survival_ratios(model_f3, data_all, data_f3, [6.0])
# fit_f44 = fit_survival_ratios(model_f44, data_all, data_f44, [0.8, 0.6])
const fit_loss = fit_survival(loss_model, data_all, [0.9, 2500])

@show Unc.(fit_f3.param, fit_f3.unc)
# @show fit_f44.param
@show Unc.(fit_loss.param, fit_loss.unc)

figure()
NaCsPlot.plot_survival_data(data_f3, fmt="C0o-", label="F=3")
NaCsPlot.plot_survival_data(data_f44, fmt="C1o-", label="4, 4")
NaCsPlot.plot_survival_data(data_all, fmt="C2o", label="Total")
plot(fit_loss.plotx, fit_loss.ploty, "C2")
grid()
xlim([0, 2500])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_raw")

figure()
plot_survival_ratios(data_all, data_f3, fmt="C0o", label="F=3")
plot(fit_f3.plotx, fit_f3.ploty, "C0")
plot_survival_ratios(data_all, data_f44, fmt="C1o-", label="4, 4")
# plot(fit_f44.plotx, fit_f44.ploty, "C1")
grid()
xlim([0, 2500])
ylim([0, 1])
legend()
title("Normalized")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ratio")

scatter_init = rates_trap * init_mf * fit_f3.param[1] * 1000
@show scatter_init sum(scatter_init)
Ω_trap = Unc.(√(fit_f3.param[1]), fit_f3.unc[1] / √(fit_f3.param[1]) / 2) * √(1000) * 1e10
@show Ω_trap
@show Ω_trap / 2π

figure()
title("Simulated")
plot_ts = linspace(0, 2500, 1000)
function compute_mfs(ts)
    A_trap = rates_to_A(rates_trap)
    res = Matrix{Float64}(16, length(ts))
    for i in 1:length(ts)
        t = ts[i]
        res[:, i] = exp(A_trap * (t * fit_f3.param[1])) * init_mf # expm
    end
    return res
end
mfs = compute_mfs(plot_ts)
plot(plot_ts, mfs[9, :], "C0", label="4, 4")
plot(plot_ts, mfs[8, :], "C2", label="4, 3")
plot(plot_ts, mfs[16, :], "C1", label="3, 3")
plot(plot_ts, mfs[15, :], "C3", label="3, 2")
grid()
xlim([0, 2500])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_sim")

figure()
title("Simulated with loss")
plot_ts = linspace(0, 2500, 1000)
function compute_mfs_loss(ts, loss)
    A_trap = rates_to_A(rates_trap) * fit_f3.param[1] - diagm(0=>loss) - I / fit_loss.param[2]
    res = Matrix{Float64}(16, length(ts))
    for i in 1:length(ts)
        t = ts[i]
        res[:, i] = exp(A_trap * t) * init_mf # expm
    end
    return res
end
mfs_loss = compute_mfs_loss(plot_ts, [1, 1, 1, 1, 1, 1, 1, 1, 0.0,
                                      1, 1, 1, 1, 1, 1, 0.0] * 1000)
plot(plot_ts, mfs_loss[9, :], "C0", label="4, 4")
# plot(plot_ts, mfs_loss[8, :], "C2", label="4, 3")
plot(plot_ts, mfs_loss[16, :], "C1", label="3, 3")
# plot(plot_ts, mfs_loss[15, :], "C3", label="3, 2")
plot(plot_ts, sum(mfs_loss, 1)[1, :], "C4", label="Total")
plot(plot_ts, exp.(-plot_ts ./ fit_loss.param[2]), "C5", label="Single body")
grid()
xlim([0, 2500])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_sim_loss")


NaCsPlot.maybe_show()
