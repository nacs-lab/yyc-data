#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180909_145616.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, -2), (3, -4)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((-1, 2), (-3, 4)))
data_nacs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2,), (3, 4,)))

const spec_a = OrderedDict(
    :lifetime=>[0, 10, 20, 50, 100, 200, 400] ./ 1000,
)

const split_na_a = NaCsData.split_data(data_na_a, spec_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)
const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)

data_na_lifetime = split_na_a[:lifetime]
data_cs_lifetime = split_cs_a[:lifetime]
data_nacs_lifetime = split_nacs_a[:lifetime]

const prefix = joinpath(@__DIR__, "imgs", "data_20180909_145616_2body_lifetime")

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
lifetime_model2(x, p) = p[1] .* exp.(.-x ./ p[2]) .+ p[3]

const plotx = linspace(0, 0.5, 1001)
fit_na = fit_survival(lifetime_model, data_na_lifetime, [0.95, 1.0], plotx=plotx)
fit_cs = fit_survival(lifetime_model, data_cs_lifetime, [0.95, 1.0], plotx=plotx)
fit_nacs = fit_survival(lifetime_model, data_nacs_lifetime, [0.95, 0.2], plotx=plotx)
fit_nacs2 = fit_survival(lifetime_model2, data_nacs_lifetime, [0.95, 0.2, 0.1], plotx=plotx)

function fit_text2(xy1, xy2, name, fit; kws...)
    t = Unc(fit.param[2], fit.unc[2])
    offset = Unc(fit.param[3], fit.unc[3])
    text(xy1..., "\$\\tau_{$name}=$(t)s\$"; kws...)
    text(xy2..., "\$off=$(offset)\$"; kws...)
end

figure()
NaCsPlot.plot_survival_data(data_cs_lifetime, fmt="C0.")
plot(fit_cs.plotx, fit_cs.ploty, "C0-", label="Cs")
NaCsPlot.plot_survival_data(data_na_lifetime, fmt="C1.")
plot(fit_na.plotx, fit_na.ploty, "C1-", label="Na")
NaCsPlot.plot_survival_data(data_nacs_lifetime, fmt="C2.")
plot(fit_nacs.plotx, fit_nacs.ploty, "C2-", label="Na+Cs")
plot(fit_nacs2.plotx, fit_nacs2.ploty, "C3-", label="With offset")
text(0.3, 0.9, "\$\\tau_{Cs}=$(Unc(fit_cs.param[2], fit_cs.unc[2]))s\$", color="C0")
text(0.05, 0.7, "\$\\tau_{Na}=$(Unc(fit_na.param[2], fit_na.unc[2]))s\$", color="C1")
text(0.03, 0.05, "\$\\tau_{Na+Cs}=$(Unc(fit_nacs.param[2], fit_nacs.unc[2]))s\$", color="C2")
fit_text2((0.32, 0.25), (0.32, 0.17), "", fit_nacs2, color="C3")
grid()
legend(fontsize=16)
xlim([0, 0.55])
ylim([0, 1])
title("Lifetime")
xlabel("Time (s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
