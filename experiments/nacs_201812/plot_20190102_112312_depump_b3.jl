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

const inames = ["data_20190102_112312.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_cs = [OrderedDict(:align=>([0, 0.5, 1, 2, 5, 10, 20, 40],
                                       [0, 0.5, 1, 2, 5, 10, 20, 40],
                                       [0, 0.5, 1, 2, 5, 10, 20, 40]),
                              :mis=>([0, 0.02, 0.05, 0.1, 0.2, 0.5],
                                     [0, 0.02, 0.05, 0.1, 0.2, 0.5],
                                     [0, 0.02, 0.05, 0.1, 0.2, 0.5]))]
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

const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_cs_align = datas_cs[1][:align]
data_cs_mis = datas_cs[1][:mis]

fit_cs_align = [fit_survival(model_exp, data, [-0.9, 5, 0.9]) for data in data_cs_align]
fit_cs_mis = [fit_survival(model_exp, data, [-0.9, 0.1, 0.9]) for data in data_cs_mis]

const prefix = joinpath(@__DIR__, "imgs", "data_20190102_112312_depump")

for i in 1:3
    figure()
    NaCsPlot.plot_survival_data(data_cs_align[i], fmt="C0.", label="\$0\\!^\\circ{}\$")
    plot(fit_cs_align[i].plotx, fit_cs_align[i].ploty, "C0")
    NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v * 100, data_cs_mis[i]),
                                fmt="C1.", label="\$45\\!^\\circ{}_{t\\times100}\$")
    plot(fit_cs_mis[i].plotx .* 100, fit_cs_mis[i].ploty, "C1")
    text(9, 0.42, "\$\\tau_{0\\!^\\circ{}}=$(fit_cs_align[i].uncs[2])ms\$", color="C0")
    text(9, 0.31, "\$\\tau_{45\\!^\\circ{}}=$(fit_cs_mis[i].uncs[2])ms\$", color="C1")
    text(9, 0.18, "\$r=$(fit_cs_align[i].uncs[2] / fit_cs_mis[i].uncs[2])\$")
    legend(fontsize="small")
    grid()
    ylim([0, 1])
    xlim([0, 55])
    title("Cs depumping \$B_$i\$")
    xlabel("Time (ms)")
    ylabel("Survival")
    NaCsPlot.maybe_save("$(prefix)_cs_b$i")
end

figure()
for i in 1:3
    NaCsPlot.plot_survival_data(data_cs_align[i], fmt="C$(i - 1).", label="\$B_$i\$")
    plot(fit_cs_align[i].plotx, fit_cs_align[i].ploty, "C$(i - 1)")
end
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 55])
title("Cs depumping \$0^\\circ\$")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_align")

figure()
for i in 1:3
    NaCsPlot.plot_survival_data(data_cs_mis[i], fmt="C$(i - 1).", label="\$B_$i\$")
    plot(fit_cs_mis[i].plotx, fit_cs_mis[i].ploty, "C$(i - 1)")
end
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 0.55])
title("Cs depumping \$45^\\circ\$")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_mis")

NaCsPlot.maybe_show()
