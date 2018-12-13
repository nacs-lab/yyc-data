#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181211_192658.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [0:8.0:160]
const specs_cs = [0:40.0:800]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function fit_survival(model, data, p0; plotx=nothing, plot_lo=nothing, plot_hi=nothing,
                      use_unc=true, plot_scale=1.1)
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
        if plot_lo === nothing
            plot_lo = mid - span * plot_scale / 2
            if plot_lo * lo <= 0
                plot_lo = 0
            end
        end
        if plot_hi === nothing
            plot_hi = mid + span * plot_scale / 2
            if plot_hi * hi <= 0
                plot_hi = 0
            end
        end
        plotx = linspace(plot_lo, plot_hi, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios[:, 2], uncs[:, 2].^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_na = datas_na[1]
data_cs = datas_cs[1]

model_cos1(x, p) = p[1] .+ p[2] .* (1 .+ cos.(π .* x ./ p[3]))
model_cos2(x, p) = p[1] .+ p[2] .* (1 .+ cos.(π .* x ./ p[3])) .+
    p[4] .* (1 .+ cos.(π .* x ./ p[5]) .* exp.(.-x ./ p[6]))

fit_na = fit_survival(model_cos1, data_na, [0.1, 0.7, 60], use_unc=false)
# @show fit_na.uncs
fit_cs = fit_survival(model_cos2, data_cs, [0.1, 0.5, 400, 0.2, 100, 200], use_unc=false)
# @show fit_cs.uncs

const prefix = joinpath(@__DIR__, "imgs", "data_20181211_192658_reverse_heating_rabi")

figure()
NaCsPlot.plot_survival_data(data_na, fmt="C0.")
plot(fit_na.plotx, fit_na.ploty, "C0-")
grid()
ylim([0, 1])
xlim([0, 170])
text(25, 0.65, "\$t_\\pi=$(fit_na.uncs[3])\$us", color="C1")
title("Na heating")
xlabel("Raman time (us)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

figure()
NaCsPlot.plot_survival_data(data_cs, fmt="C0.")
plot(fit_cs.plotx, fit_cs.ploty, "C0-")
grid()
ylim([0, 1])
xlim([0, 900])
text(70, 0.6, "\$t_{\\pi1}=$(fit_cs.uncs[3])\$us\n" *
     "\$t_{\\pi2}=$(fit_cs.uncs[5])\$us\n" *
     "\$\\tau_{2}=$(fit_cs.uncs[6])\$us", color="C1")
title("Cs heating")
xlabel("Raman time (us)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
