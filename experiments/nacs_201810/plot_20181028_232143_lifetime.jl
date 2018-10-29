#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181028_232143.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_a = OrderedDict(
    :nocool=>[0, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
    :cool=>[0, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
)

split_na = NaCsData.split_data(data_na_a, spec_a)
split_cs = NaCsData.split_data(data_cs_a, spec_a)

data_na_nocool = split_na[:nocool]
data_na_cool = split_na[:cool]
data_cs_nocool = split_cs[:nocool]
data_cs_cool = split_cs[:cool]

const prefix = joinpath(@__DIR__, "imgs", "data_20181028_232143_lifetime")

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

lifetime_model(x, p) = p[1] .* exp.(.-x ./ p[2])
fit_na_nocool = fit_survival(lifetime_model, data_na_nocool, [0.95, 1.4],
                             plotx=linspace(0, 2.5, 1051))
fit_na_cool = fit_survival(lifetime_model, data_na_cool, [0.95, 1.4],
                           plotx=linspace(0, 2.5, 1051))
fit_cs_nocool = fit_survival(lifetime_model, data_cs_nocool, [0.95, 4.4],
                             plotx=linspace(0, 10.5, 1051))
fit_cs_cool = fit_survival(lifetime_model, data_cs_cool, [0.95, 4.4],
                           plotx=linspace(0, 10.5, 1051))

figure()
NaCsPlot.plot_survival_data(data_na_nocool, fmt="C0.")
plot(fit_na_nocool.plotx, fit_na_nocool.ploty, "C0-", label="No cool")
NaCsPlot.plot_survival_data(data_na_cool, fmt="C1.")
plot(fit_na_cool.plotx, fit_na_cool.ploty, "C1-", label="With cool")
text(1, 0.55, "\$\\tau_{nocool}=$(Unc(fit_na_nocool.param[2], fit_na_nocool.unc[2]))s\$",
     color="C0")
text(1, 0.45, "\$\\tau_{cool}=$(Unc(fit_na_cool.param[2], fit_na_cool.unc[2]))s\$",
     color="C1")
grid()
legend()
ylim([0, 1])
xlim([0, 2.5])
title("Na lifetime")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

figure()
NaCsPlot.plot_survival_data(data_cs_nocool, fmt="C0.")
plot(fit_cs_nocool.plotx, fit_cs_nocool.ploty, "C0-", label="No cool")
NaCsPlot.plot_survival_data(data_cs_cool, fmt="C1.")
plot(fit_cs_cool.plotx, fit_cs_cool.ploty, "C1-", label="With cool")
text(5, 0.55, "\$\\tau_{nocool}=$(Unc(fit_cs_nocool.param[2], fit_cs_nocool.unc[2]))s\$",
     color="C0")
text(5, 0.45, "\$\\tau_{cool}=$(Unc(fit_cs_cool.param[2], fit_cs_cool.unc[2]))s\$",
     color="C1")
grid()
legend()
ylim([0, 1])
title("Cs lifetime")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
