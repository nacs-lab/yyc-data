#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190101_160615.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [0:2:40]
const specs_cs = [0:2:40]
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

model_sin0(x, p) = p[1] .+ p[2] .* sin.(2π .* x .* p[3] .+ π / 2)

data_na = datas_na[1]
data_cs = datas_cs[1]

fit_na = fit_survival(model_sin0, data_na, [0.45, 0.45, 0.036])
fit_cs = fit_survival(model_sin0, data_cs, [0.45, 0.45, 0.02])

const prefix = joinpath(@__DIR__, "imgs", "data_20190101_160615_coprop_back_time")

figure()
NaCsPlot.plot_survival_data(data_na, fmt="C0.")
plot(fit_na.plotx, fit_na.ploty)
text(3, 0.88, "\$t_\\pi=$(1 / 2 / fit_na.uncs[3])\\mu s\$", color="C0")
grid()
ylim([0, 1])
xlim([0, 45])
title("Na coprop")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

figure()
NaCsPlot.plot_survival_data(data_cs, fmt="C0.")
plot(fit_cs.plotx, fit_cs.ploty)
text(10, 0.85, "\$t_\\pi=$(1 / 2 / fit_cs.uncs[3])\\mu s\$", color="C0")
grid()
ylim([0, 1])
xlim([0, 45])
title("Cs coprop")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
