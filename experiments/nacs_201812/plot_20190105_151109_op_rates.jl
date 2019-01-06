#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
using NaCsCalc.Atomic: all_scatter_D



const inames = ["data_20190105_151109.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [([0.0, 2, 4, 10, 20, 40, 100],
                   [0.0, 2, 4, 10, 20, 40, 100])]
const specs_cs = [([0.0, 2, 4, 10, 20, 40, 100],
                   [0.0, 2, 4, 10, 20, 40, 100])]
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

model_exp(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_cs_43 = datas_cs[1][1]
data_cs_33 = datas_cs[1][2]
data_na_21 = datas_na[1][1]
data_na_11 = datas_na[1][2]

fit_cs_43 = fit_survival(model_exp, data_cs_43, [-0.5, 10, 0.5])
fit_cs_33 = fit_survival(model_exp, data_cs_33, [0.9, 10, 0])
fit_na_21 = fit_survival(model_exp, data_na_21, [-0.5, 10, 0.5])
fit_na_11 = fit_survival(model_exp, data_na_11, [0.9, 10, 0])

const prefix = joinpath(@__DIR__, "imgs", "data_20190105_151109_op_rates")

figure()
NaCsPlot.plot_survival_data(data_cs_43, fmt="C0.", label="4, 4")
plot(fit_cs_43.plotx, fit_cs_43.ploty, "C0")
NaCsPlot.plot_survival_data(data_cs_33, fmt="C1.", label="4, 3")
plot(fit_cs_33.plotx, fit_cs_33.ploty, "C1")
text(10, 0.62, "\$\\tau_{4, 4}=$(fit_cs_43.uncs[2])\\mu s\$", color="C0")
text(10, 0.51, "\$\\tau_{4, 3}=$(fit_cs_33.uncs[2])\\mu s\$", color="C1")
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 30])
title("Cs OP rates")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

figure()
NaCsPlot.plot_survival_data(data_na_21, fmt="C0.", label="2, 2")
plot(fit_na_21.plotx, fit_na_21.ploty, "C0")
NaCsPlot.plot_survival_data(data_na_11, fmt="C1.", label="2, 1")
plot(fit_na_11.plotx, fit_na_11.ploty, "C1")
text(10, 0.42, "\$\\tau_{2, 2}=$(fit_na_21.uncs[2])\\mu s\$", color="C0")
text(10, 0.31, "\$\\tau_{2, 1}=$(fit_na_11.uncs[2])\\mu s\$", color="C1")
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 30])
title("Na OP rates")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

NaCsPlot.maybe_show()
