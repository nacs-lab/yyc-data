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

const iname_a = joinpath(@__DIR__, "data", "data_20180208_194456.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function gen_selector(na)
    function selector(logicals)
        @assert size(logicals, 2) == 1
        if logicals[2, 1] == 0 || logicals[1, 1] != na
            return Int[1, 0, 0]
        end
        return Int[1, 1, logicals[4, 1]]
    end
end

data_single_a = NaCsData.select_count(params_a, logicals_a, gen_selector(false))
data_both_a = NaCsData.select_count(params_a, logicals_a, gen_selector(true))

const spec_a = OrderedDict(
    :total=>[0, 12.5, 25, 50, 100, 200, 400, 800, 1600],
    :f3=>[0, 12.5, 25, 50, 100, 200, 400, 800, 1600],
    :total_rb=>[0, 12.5, 25, 50, 100, 200, 400, 800, 1600],
    :f3_rb=>[0, 12.5, 25, 50, 100, 200, 400, 800, 1600],
    :total_pa=>[0, 12.5, 25, 50, 100, 200, 400, 800, 1600] ./ 16,
    :f3_pa=>[0, 12.5, 25, 50, 100, 200, 400, 800, 1600] ./ 16,
    :total_rb_pa=>[0, 12.5, 25, 50, 100, 200, 400, 800, 1600] ./ 16,
    :f3_rb_pa=>[0, 12.5, 25, 50, 100, 200, 400, 800, 1600] ./ 16,
)

const split_single_a = NaCsData.split_data(data_single_a, spec_a)
const split_both_a = NaCsData.split_data(data_both_a, spec_a)

const data_single_total = split_single_a[:total]
const data_single_f3 = split_single_a[:f3]
const data_single_total_rb = split_single_a[:total_rb]
const data_single_f3_rb = split_single_a[:f3_rb]
const data_single_total_pa = split_single_a[:total_pa]
const data_single_f3_pa = split_single_a[:f3_pa]
const data_single_total_rb_pa = split_single_a[:total_rb_pa]
const data_single_f3_rb_pa = split_single_a[:f3_rb_pa]

const data_both_total = split_both_a[:total]
const data_both_f3 = split_both_a[:f3]
const data_both_total_rb = split_both_a[:total_rb]
const data_both_f3_rb = split_both_a[:f3_rb]
const data_both_total_pa = split_both_a[:total_pa]
const data_both_f3_pa = split_both_a[:f3_pa]
const data_both_total_rb_pa = split_both_a[:total_rb_pa]
const data_both_f3_rb_pa = split_both_a[:f3_rb_pa]

const prefix = joinpath(@__DIR__, "imgs", "data_20180208_194456_trap_pa_scatter")

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
const rates_trap_rb = all_scatters(cs_atom, 1e10, 299792458 / 976e-9, (0.25, 0.5, 0.25))
const rates_pa = all_scatters(cs_atom, 1e9, 351553e9, (0.4, 0.3, 0.3))
const rates_pa_rb = all_scatters(cs_atom, 1e9, 351553e9, (0.8, 0, 0.2))

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
    (x, p) -> propagate_f.(x .* p[1])
end

function gen_model2(rates1, rates2, init, idxs)
    A1 = rates_to_A(rates1)
    A2 = rates_to_A(rates2)
    (x, p) -> propagate.((A1 .+ A2 .* p[1],), (init,), x, (idxs,))
end

function loss_model(x, p)
    return p[1] .* exp.(-x ./ p[2])
end

const init_mf = [0, 0, 0, 0, 0, 0, 0, 0, 1.0,
                 0, 0, 0, 0, 0, 0, 0]
const model_f3 = gen_model(rates_trap, init_mf, (10, 11, 12, 13, 14, 15, 16))
const model_f3_rb = gen_model(rates_trap_rb, init_mf, (10, 11, 12, 13, 14, 15, 16))

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

const fit_single_total = fit_survival(loss_model, data_single_total, [0.9, 2500])
const fit_single_total_rb = fit_survival(loss_model, data_single_total_rb, [0.9, 2500])
const fit_single_total_pa = fit_survival(loss_model, data_single_total_pa, [0.9, 2500])
const fit_single_total_rb_pa = fit_survival(loss_model, data_single_total_rb_pa, [0.9, 2500])

const fit_f3 = fit_survival_ratios(model_f3, data_single_total, data_single_f3, [6.0])
const fit_f3_rb = fit_survival_ratios(model_f3_rb, data_single_total_rb, data_single_f3_rb, [6.0])

const model_f3_pa = gen_model2(rates_trap * fit_f3.param[1], rates_pa, init_mf,
                               (10, 11, 12, 13, 14, 15, 16))
const model_f3_rb_pa = gen_model2(rates_trap_rb * fit_f3.param[1], rates_pa_rb,
                                  init_mf, (10, 11, 12, 13, 14, 15, 16))

const fit_f3_pa = fit_survival_ratios(model_f3_pa, data_single_total_pa, data_single_f3_pa,
                                      [1.0])
const fit_f3_rb_pa = fit_survival_ratios(model_f3_rb_pa, data_single_total_rb_pa,
                                         data_single_f3_rb_pa, [1.0])

@show Unc.(fit_single_total.param, fit_single_total.unc)
@show Unc.(fit_single_total_rb.param, fit_single_total_rb.unc)
@show Unc.(fit_single_total_pa.param, fit_single_total_pa.unc)
@show Unc.(fit_single_total_rb_pa.param, fit_single_total_rb_pa.unc)

@show Unc.(fit_f3.param, fit_f3.unc)
@show Unc.(fit_f3_rb.param, fit_f3_rb.unc)
@show Unc.(fit_f3_pa.param, fit_f3_pa.unc)
@show Unc.(fit_f3_rb_pa.param, fit_f3_rb_pa.unc)

figure()
title("Normal B field")
NaCsPlot.plot_survival_data(data_single_total, fmt="C0o", label="Total")
plot(fit_single_total.plotx, fit_single_total.ploty, "C0")
NaCsPlot.plot_survival_data(data_single_f3, fmt="C1o-", label="F=3")
NaCsPlot.plot_survival_data(data_both_total, fmt="C2o-", label="Both")
NaCsPlot.plot_survival_data(data_both_f3, fmt="C3o-", label="Both F=3")
grid()
xlim([0, 2000])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_raw")

figure()
title("Rotated B field")
NaCsPlot.plot_survival_data(data_single_total_rb, fmt="C0o", label="Total")
plot(fit_single_total_rb.plotx, fit_single_total_rb.ploty, "C0")
NaCsPlot.plot_survival_data(data_single_f3_rb, fmt="C1o-", label="F=3")
NaCsPlot.plot_survival_data(data_both_total_rb, fmt="C2o-", label="Both")
NaCsPlot.plot_survival_data(data_both_f3_rb, fmt="C3o-", label="Both F=3")
grid()
xlim([0, 2000])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_raw_rb")

figure()
title("Normal B field with PA")
NaCsPlot.plot_survival_data(data_single_total_pa, fmt="C0o", label="Total")
plot(fit_single_total_pa.plotx, fit_single_total_pa.ploty, "C0")
NaCsPlot.plot_survival_data(data_single_f3_pa, fmt="C1o-", label="F=3")
NaCsPlot.plot_survival_data(data_both_total_pa, fmt="C2o-", label="Both")
NaCsPlot.plot_survival_data(data_both_f3_pa, fmt="C3o-", label="Both F=3")
grid()
xlim([0, 120])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_raw_pa")

figure()
title("Rotated B field with PA")
NaCsPlot.plot_survival_data(data_single_total_rb_pa, fmt="C0o", label="Total")
plot(fit_single_total_rb_pa.plotx, fit_single_total_rb_pa.ploty, "C0")
NaCsPlot.plot_survival_data(data_single_f3_rb_pa, fmt="C1o-", label="F=3")
NaCsPlot.plot_survival_data(data_both_total_rb_pa, fmt="C2o-", label="Both")
NaCsPlot.plot_survival_data(data_both_f3_rb_pa, fmt="C3o-", label="Both F=3")
grid()
xlim([0, 120])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_raw_rb_pa")

figure()
plot_survival_ratios(data_single_total, data_single_f3, fmt="C0o", label="Normal")
plot(fit_f3.plotx, fit_f3.ploty, "C0")
plot_survival_ratios(data_single_total_rb, data_single_f3_rb, fmt="C1o", label="Rotated B")
plot(fit_f3_rb.plotx, fit_f3_rb.ploty, "C1")
plot_survival_ratios(data_both_total, data_both_f3, fmt="C2", label="Both")
plot_survival_ratios(data_both_total_rb, data_both_f3_rb, fmt="C3", label="Both Rotate B")
grid()
xlim([0, 2000])
ylim([0, 1])
legend()
title("Normalized")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ratio")

figure()
plot_survival_ratios(data_single_total_pa, data_single_f3_pa, fmt="C0o", label="Normal")
plot(fit_f3_pa.plotx, fit_f3_pa.ploty, "C0")
plot_survival_ratios(data_single_total_rb_pa, data_single_f3_rb_pa, fmt="C1o", label="Rotated B")
plot(fit_f3_rb_pa.plotx, fit_f3_rb_pa.ploty, "C1")
grid()
xlim([0, 120])
ylim([0, 1])
legend()
title("Normalized with PA")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ratio_pa")

# scatter_init = rates_trap * init_mf * fit_f3.param[1] * 1000
# @show scatter_init sum(scatter_init)
# Ω_trap = Unc.(√(fit_f3.param[1]), fit_f3.unc[1] / √(fit_f3.param[1]) / 2) * √(1000) * 1e10
# @show Ω_trap
# @show Ω_trap / 2π

# figure()
# title("Simulated")
# plot_ts = linspace(0, 2500, 1000)
# function compute_mfs(ts)
#     A_trap = rates_to_A(rates_trap)
#     res = Matrix{Float64}(16, length(ts))
#     for i in 1:length(ts)
#         t = ts[i]
#         res[:, i] = exp(A_trap * (t * fit_f3.param[1])) * init_mf # expm
#     end
#     return res
# end
# mfs = compute_mfs(plot_ts)
# plot(plot_ts, mfs[9, :], "C0", label="4, 4")
# plot(plot_ts, mfs[8, :], "C2", label="4, 3")
# plot(plot_ts, mfs[16, :], "C1", label="3, 3")
# plot(plot_ts, mfs[15, :], "C3", label="3, 2")
# grid()
# xlim([0, 2500])
# ylim([0, 1])
# legend()
# xlabel("Time (\$ms\$)")
# ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)_sim")

# figure()
# title("Simulated with loss")
# plot_ts = linspace(0, 2500, 1000)
# function compute_mfs_loss(ts, loss)
#     A_trap = rates_to_A(rates_trap) * fit_f3.param[1] - diagm(0=>loss) - I / fit_loss.param[2]
#     res = Matrix{Float64}(16, length(ts))
#     for i in 1:length(ts)
#         t = ts[i]
#         res[:, i] = exp(A_trap * t) * init_mf # expm
#     end
#     return res
# end
# mfs_loss = compute_mfs_loss(plot_ts, [1, 1, 1, 1, 1, 1, 1, 1, 0.0,
#                                       1, 1, 1, 1, 1, 1, 0.0] * 1000)
# plot(plot_ts, mfs_loss[9, :], "C0", label="4, 4")
# # plot(plot_ts, mfs_loss[8, :], "C2", label="4, 3")
# plot(plot_ts, mfs_loss[16, :], "C1", label="3, 3")
# # plot(plot_ts, mfs_loss[15, :], "C3", label="3, 2")
# plot(plot_ts, sum(mfs_loss, 1)[1, :], "C4", label="Total")
# plot(plot_ts, exp.(-plot_ts ./ fit_loss.param[2]), "C5", label="Single body")
# grid()
# xlim([0, 2500])
# ylim([0, 1])
# legend()
# xlabel("Time (\$ms\$)")
# ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)_sim_loss")

NaCsPlot.maybe_show()
