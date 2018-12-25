#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181220_154644.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [i->i <= 8670 || i > 9180]
const specs = [0:8.0:400]
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

model_sin(x, p) = p[1] .+ p[2] .* sin.(2π .* x .* p[3] .+ p[4])

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs)

data_na = datas_na[1]
data_cs = datas_cs[1]

fit_na = fit_survival(model_sin, data_na, [0.5, 0.12, 0.010, 1])
fit_cs = fit_survival(model_sin, data_cs, [0.5, 0.3, 0.006, -2])

const prefix = joinpath(@__DIR__, "imgs", "data_20181220_154644_diff_ramsey")

figure()
NaCsPlot.plot_survival_data(data_na, fmt="C0.", label="Na")
plot(fit_na.plotx, fit_na.ploty, "C0-")
NaCsPlot.plot_survival_data(data_cs, fmt="C1.", label="Cs")
plot(fit_cs.plotx, fit_cs.ploty, "C1-")
text(30, 0.85, "\$\\nu_{Na}=$(fit_na.uncs[3] * 1000)\$kHz", color="C0")
text(90, 0.04, "\$\\nu_{Cs}=$(fit_cs.uncs[3] * 1000)\$kHz", color="C1")
legend()
grid()
ylim([0, 1])
title("X-Y diff Ramsey")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
