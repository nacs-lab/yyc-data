#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180831_175841.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_a = OrderedDict(
    :lifetime=>[0, 0.2, 0.5, 1, 2, 5, 10],
    :loadp=>4.0:20.0,
)

const split_na_a = NaCsData.split_data(data_na_a, spec_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)

data_na_lifetime = split_na_a[:lifetime]
data_cs_lifetime = split_cs_a[:lifetime]
data_cs_loadp = split_cs_a[:loadp]

const prefix = joinpath(@__DIR__, "imgs", "data_20180831_175841_csload")

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

fit_na = fit_survival(lifetime_model, data_na_lifetime, [0.95, 4.4], plotx=linspace(0, 11, 1001))
fit_cs = fit_survival(lifetime_model, data_cs_lifetime, [0.95, 4.4], plotx=linspace(0, 11, 1001))

figure()
NaCsPlot.plot_survival_data(data_na_lifetime, fmt="C0.")
plot(fit_na.plotx, fit_na.ploty, "C0-", label="Na")
NaCsPlot.plot_survival_data(data_cs_lifetime, fmt="C1.")
plot(fit_cs.plotx, fit_cs.ploty, "C1-", label="Cs")
text(5, 0.55, "\$\\tau_{Na}=$(Unc(fit_na.param[2], fit_na.unc[2]))s\$", color="C0")
text(5, 0.45, "\$\\tau_{Cs}=$(Unc(fit_cs.param[2], fit_cs.unc[2]))s\$", color="C1")
grid()
legend()
xlim([0, 11.5])
ylim([0, 1])
title("Lifetime")
xlabel("Time (s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_lifetime")

figure()
NaCsPlot.plot_loading_data(data_cs_loadp, fmt="C0.-", label="Loading")
NaCsPlot.plot_survival_data(data_cs_loadp, fmt="C1.-", label="Surviving")
grid()
legend()
ylim([0, 1])
title("Cs loading")
xlabel("Tweezer power (mW)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_loadp")

NaCsPlot.maybe_show()
