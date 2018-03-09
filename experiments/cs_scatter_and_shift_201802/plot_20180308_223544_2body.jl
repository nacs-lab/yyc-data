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

const iname_a = joinpath(@__DIR__, "data", "data_20180308_223544.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[2, 1] == 0 || logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[4, 1] != 0 && logicals[3, 1] != 0]
end

function gen_single_selector(na)
    function selector(logicals)
        @assert size(logicals, 2) == 1
        if logicals[2, 1] != !na || logicals[1, 1] != na
            return Int[1, 0, 0]
        end
        return Int[1, 1, logicals[4, 1] == !na && logicals[3, 1] == na]
    end
end

data_both_a = NaCsData.select_count(params_a, logicals_a, selector)
data_na_a = NaCsData.select_count(params_a, logicals_a, gen_single_selector(true))
data_cs_a = NaCsData.select_count(params_a, logicals_a, gen_single_selector(false))

const spec_a = OrderedDict(
    :op_14=>[0, 20, 50, 100, 200, 400],
    :noop_14=>[0, 20, 50, 100, 200, 400],
    :op_29=>[0, 20, 50, 100, 200, 400],
    :noop_29=>[0, 20, 50, 100, 200, 400],
    :op_7=>[0, 20, 50, 100, 200, 400],
    :noop_7=>[0, 20, 50, 100, 200, 400],
)

const split_both_a = NaCsData.split_data(data_both_a, spec_a)
const split_na_a = NaCsData.split_data(data_na_a, spec_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)

const data_both_op_14 = split_both_a[:op_14]
const data_both_noop_14 = split_both_a[:noop_14]
const data_both_op_29 = split_both_a[:op_29]
const data_both_noop_29 = split_both_a[:noop_29]
const data_both_op_7 = split_both_a[:op_7]
const data_both_noop_7 = split_both_a[:noop_7]

const data_na_op_14 = split_na_a[:op_14]
const data_na_noop_14 = split_na_a[:noop_14]
const data_na_op_29 = split_na_a[:op_29]
const data_na_noop_29 = split_na_a[:noop_29]
const data_na_op_7 = split_na_a[:op_7]
const data_na_noop_7 = split_na_a[:noop_7]

const data_cs_op_14 = split_cs_a[:op_14]
const data_cs_noop_14 = split_cs_a[:noop_14]
const data_cs_op_29 = split_cs_a[:op_29]
const data_cs_noop_29 = split_cs_a[:noop_29]
const data_cs_op_7 = split_cs_a[:op_7]
const data_cs_noop_7 = split_cs_a[:noop_7]

const prefix = joinpath(@__DIR__, "imgs", "data_20180308_223544_2body")

function loss_model(x, p)
    return p[1] .* exp.(-x .* p[2])
end

# function get_survival_data(data, z=1.0)
#     params, ratios, uncs = NaCsData.get_values(data, z)
#     perm = sortperm(params)
#     params = params[perm]
#     return params, ratios[perm, 2], uncs[perm, 2]
# end

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

const fit_noop_7 = fit_survival(loss_model, data_both_noop_7, [1.0, 0.001])
const fit_noop_14 = fit_survival(loss_model, data_both_noop_14, [1.0, 0.002])
const fit_noop_29 = fit_survival(loss_model, data_both_noop_29, [1.0, 0.005])
const fit_op_7 = fit_survival(loss_model, data_both_op_7, [1.0, 0.001])
const fit_op_14 = fit_survival(loss_model, data_both_op_14, [1.0, 0.002])
const fit_op_29 = fit_survival(loss_model, data_both_op_29, [1.0, 0.005])

const fit_na_noop_7 = fit_survival(loss_model, data_na_noop_7, [1.0, 0.001])
const fit_na_noop_14 = fit_survival(loss_model, data_na_noop_14, [1.0, 0.001])
const fit_na_noop_29 = fit_survival(loss_model, data_na_noop_29, [1.0, 0.001])
const fit_na_op_7 = fit_survival(loss_model, data_na_op_7, [1.0, 0.001])
const fit_na_op_14 = fit_survival(loss_model, data_na_op_14, [1.0, 0.001])
const fit_na_op_29 = fit_survival(loss_model, data_na_op_29, [1.0, 0.001])

const fit_cs_noop_7 = fit_survival(loss_model, data_cs_noop_7, [1.0, 0.001])
const fit_cs_noop_14 = fit_survival(loss_model, data_cs_noop_14, [1.0, 0.001])
const fit_cs_noop_29 = fit_survival(loss_model, data_cs_noop_29, [1.0, 0.001])
const fit_cs_op_7 = fit_survival(loss_model, data_cs_op_7, [1.0, 0.001])
const fit_cs_op_14 = fit_survival(loss_model, data_cs_op_14, [1.0, 0.001])
const fit_cs_op_29 = fit_survival(loss_model, data_cs_op_29, [1.0, 0.001])

# @show Unc.(fit_noop_7.param, fit_noop_7.unc)
# @show Unc.(fit_noop_14.param, fit_noop_14.unc)
# @show Unc.(fit_noop_29.param, fit_noop_29.unc)
# @show Unc.(fit_op_7.param, fit_op_7.unc)
# @show Unc.(fit_op_14.param, fit_op_14.unc)
# @show Unc.(fit_op_29.param, fit_op_29.unc)

figure()
title("No OP 1-body loss")
NaCsPlot.plot_survival_data(data_na_noop_7, fmt="C0.", label="Na 7")
plot(fit_na_noop_7.plotx, fit_na_noop_7.ploty, "C0")
NaCsPlot.plot_survival_data(data_na_noop_14, fmt="C1.", label="Na 14.3")
plot(fit_na_noop_14.plotx, fit_na_noop_14.ploty, "C1")
NaCsPlot.plot_survival_data(data_na_noop_29, fmt="C2.", label="Na 29")
plot(fit_na_noop_29.plotx, fit_na_noop_29.ploty, "C2")
NaCsPlot.plot_survival_data(data_cs_noop_7, fmt="C3.", label="Cs 7")
plot(fit_cs_noop_7.plotx, fit_cs_noop_7.ploty, "C3")
NaCsPlot.plot_survival_data(data_cs_noop_14, fmt="C4.", label="Cs 14.3")
plot(fit_cs_noop_14.plotx, fit_cs_noop_14.ploty, "C4")
NaCsPlot.plot_survival_data(data_cs_noop_29, fmt="C5.", label="Cs 29")
plot(fit_cs_noop_29.plotx, fit_cs_noop_29.ploty, "C5")
grid()
xlim([0, 500])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_single_noop")

figure()
title("With OP 1-body loss")
NaCsPlot.plot_survival_data(data_na_op_7, fmt="C0.", label="Na 7")
plot(fit_na_op_7.plotx, fit_na_op_7.ploty, "C0")
NaCsPlot.plot_survival_data(data_na_op_14, fmt="C1.", label="Na 14.3")
plot(fit_na_op_14.plotx, fit_na_op_14.ploty, "C1")
NaCsPlot.plot_survival_data(data_na_op_29, fmt="C2.", label="Na 29")
plot(fit_na_op_29.plotx, fit_na_op_29.ploty, "C2")
NaCsPlot.plot_survival_data(data_cs_op_7, fmt="C3.", label="Cs 7")
plot(fit_cs_op_7.plotx, fit_cs_op_7.ploty, "C3")
NaCsPlot.plot_survival_data(data_cs_op_14, fmt="C4.", label="Cs 14.3")
plot(fit_cs_op_14.plotx, fit_cs_op_14.ploty, "C4")
NaCsPlot.plot_survival_data(data_cs_op_29, fmt="C5.", label="Cs 29")
plot(fit_cs_op_29.plotx, fit_cs_op_29.ploty, "C5")
grid()
xlim([0, 500])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_single_op")

figure()
title("No OP 2-body loss")
NaCsPlot.plot_survival_data(data_both_noop_7, fmt="C0.", label="7")
plot(fit_noop_7.plotx, fit_noop_7.ploty, "C0")
NaCsPlot.plot_survival_data(data_both_noop_14, fmt="C1.", label="14.3")
plot(fit_noop_14.plotx, fit_noop_14.ploty, "C1")
NaCsPlot.plot_survival_data(data_both_noop_29, fmt="C2.", label="29")
plot(fit_noop_29.plotx, fit_noop_29.ploty, "C2")
grid()
xlim([0, 500])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_noop")

figure()
title("With OP 2-body loss")
NaCsPlot.plot_survival_data(data_both_op_7, fmt="C0.", label="7")
plot(fit_op_7.plotx, fit_op_7.ploty, "C0")
NaCsPlot.plot_survival_data(data_both_op_14, fmt="C1.", label="14.3")
plot(fit_op_14.plotx, fit_op_14.ploty, "C1")
NaCsPlot.plot_survival_data(data_both_op_29, fmt="C2.", label="29")
plot(fit_op_29.plotx, fit_op_29.ploty, "C2")
grid()
xlim([0, 500])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_op")

rate_both = Unc.([fit_noop_7.param[2], fit_noop_14.param[2], fit_noop_29.param[2]],
                 [fit_noop_7.unc[2], fit_noop_14.unc[2], fit_noop_29.unc[2]])
rate_na = Unc.([fit_na_noop_7.param[2], fit_na_noop_14.param[2], fit_na_noop_29.param[2]],
               [fit_na_noop_7.unc[2], fit_na_noop_14.unc[2], fit_na_noop_29.unc[2]])
rate_cs = Unc.([fit_cs_noop_7.param[2], fit_cs_noop_14.param[2], fit_cs_noop_29.param[2]],
               [fit_cs_noop_7.unc[2], fit_cs_noop_14.unc[2], fit_cs_noop_29.unc[2]])

rate_2body = rate_both .- (rate_na .+ rate_cs)

param_2body = [r.a for r in rate_2body]
unc_2body = [r.s for r in rate_2body]

function rate_model(x, p)
    return p[1] .* x.^p[2]
end
fit_rate = curve_fit(rate_model, [7, 14.3, 29], param_2body, [2e-5, 2.0])
fit_rate = (param=fit_rate.param, unc=estimate_errors(fit_rate),
            plotx=linspace(6, 30, 1000),
            ploty=(rate_model((linspace(6, 30, 1000)), fit_rate.param)))
@show Unc.(fit_rate.param, fit_rate.unc)

figure()
errorbar([7, 14.3, 29], param_2body, unc_2body, fmt="C0.")
plot(fit_rate.plotx, fit_rate.ploty)
ax = gca()
ax[:set_xscale]("log", nonposx="clip")
ax[:set_yscale]("log", nonposy="clip")
text(6, 0.004, "\$r\\propto P^{$(Unc(fit_rate.param[2], fit_rate.unc[2]))}\$",
     size=30)
xlabel("Trap power (P/mW)")
ylabel("Loss rate (r/ms\$^{-1})\$")
grid()
NaCsPlot.maybe_save("$(prefix)_noop_rate")

NaCsPlot.maybe_show()
