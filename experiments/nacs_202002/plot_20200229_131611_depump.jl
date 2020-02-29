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

const inames = ["data_20200229_131611.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [([0, 0.5, 1, 2, 5, 10, 20] .* 50,
                   [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4] .* 0.2)]
const specs_cs = [([0, 0.5, 1, 2, 5, 10, 20] .* 10,
                   [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4] .* 1)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

fit_data(model, x, y, p0; kws...) =
    fit_data(model, x, y, nothing, p0; kws...)

function fit_data(model, params, ratios, uncs, p0;
                  plotx=nothing, plot_lo=nothing, plot_hi=nothing, plot_scale=1.1)
    use_unc = uncs !== nothing
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
        fit = curve_fit(model, params, ratios, uncs.^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios, p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

function fit_survival(model, data, p0; use_unc=true, kws...)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
        return fit_data(model, params, ratios[:, 2], uncs[:, 2], p0; kws...)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
        return fit_data(model, params, ratios[:, 2], p0; kws...)
    end
end

model_exp(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]

function gen_exp0(p0)
    return (x, p) -> p0 .* (1 .- exp.(x ./ -p[1]))
end

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_cs_44 = datas_cs[1][1]
data_cs_43 = datas_cs[1][2]
data_na_22 = datas_na[1][1]
data_na_21 = datas_na[1][2]

fit_cs_44 = fit_survival(model_exp, data_cs_44, [-0.9, 500, 0.9])
# fit_cs_44 = fit_survival(gen_exp0(0.95), data_cs_44, [5.0])
fit_cs_43 = fit_survival(model_exp, data_cs_43, [-0.5, 1, 0.5])
fit_na_22 = fit_survival(model_exp, data_na_22, [-0.9, 500, 0.9])
fit_na_21 = fit_survival(model_exp, data_na_21, [-0.5, 0.06, 0.5])

const prefix = joinpath(@__DIR__, "imgs", "data_20200229_131611_depump")

figure()
NaCsPlot.plot_survival_data(data_na_22, fmt="C0.", label="2, 2")
plot(fit_na_22.plotx, fit_na_22.ploty, "C0")
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v * 10000, data_na_21),
                            fmt="C1.", label="2, 1\${}_{t\\times10000}\$")
plot(fit_na_21.plotx .* 10000, fit_na_21.ploty, "C1")
text(150, 0.31, "\$\\tau_{2, 2}=$(fit_na_22.uncs[2])ms\$", color="C0")
text(200, 0.62, "\$\\tau_{2, 1}=$(fit_na_21.uncs[2] * 1000)\\mu s\$", color="C1")
text(50, 0.03, "\$r=$(fit_na_22.uncs[2] / fit_na_21.uncs[2])\$")
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 1100])
title("Na depumping")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

figure()
NaCsPlot.plot_survival_data(data_cs_44, fmt="C0.", label="4, 4")
plot(fit_cs_44.plotx, fit_cs_44.ploty, "C0")
NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v * 500, data_cs_43),
                            fmt="C1.", label="4, 3\${}_{t\\times500}\$")
plot(fit_cs_43.plotx .* 500, fit_cs_43.ploty, "C1")
text(80, 0.6, "\$\\tau_{4, 4}=$(fit_cs_44.uncs[2])ms\$", color="C0")
text(70, 0.35, "\$\\tau_{4, 3}=$(fit_cs_43.uncs[2] * 1000)\\mu s\$", color="C1")
text(10, 0.02, "\$r=$(fit_cs_44.uncs[2] / fit_cs_43.uncs[2])\$")
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 240])
title("Cs depumping")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
