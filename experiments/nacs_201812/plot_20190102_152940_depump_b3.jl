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

const inames = ["data_20190102_152940.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [OrderedDict(:i22=>([0, 0.5, 1, 2, 5, 10, 20, 40] ./ 2,
                                     [0, 0.5, 1, 2, 5, 10, 20, 40] ./ 2,
                                     [0, 0.5, 1, 2, 5, 10, 20, 40] ./ 2),
                              :i21=>([0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5],
                                     [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5],
                                     [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]))]
const specs_cs = [OrderedDict(:i44=>([0, 0.5, 1, 2, 5, 10, 20, 40],
                                     [0, 0.5, 1, 2, 5, 10, 20, 40],
                                     [0, 0.5, 1, 2, 5, 10, 20, 40]),
                              :i43=>([0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5],
                                     [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5],
                                     [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]))]
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
const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)

data_cs_44 = datas_cs[1][:i44]
data_cs_43 = datas_cs[1][:i43]
data_na_22 = datas_na[1][:i22]
data_na_21 = datas_na[1][:i21]

fit_cs_44 = [fit_survival(model_exp, data, [-0.9, 5, 0.9]) for data in data_cs_44]
fit_cs_43 = [fit_survival(model_exp, data, [-0.9, 0.1, 0.9]) for data in data_cs_43]
fit_na_22 = [fit_survival(model_exp, data, [-0.9, 5, 0.9]) for data in data_na_22]
fit_na_21 = [fit_survival(model_exp, data, [-0.9, 0.1, 0.9]) for data in data_na_21]

const prefix = joinpath(@__DIR__, "imgs", "data_20190102_152940_depump")

for i in 1:3
    figure()
    NaCsPlot.plot_survival_data(data_cs_44[i], fmt="C0.", label="4, 4")
    plot(fit_cs_44[i].plotx, fit_cs_44[i].ploty, "C0")
    NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v * 100, data_cs_43[i]),
                                fmt="C1.", label="4, 3\${}_{t\\times100}\$")
    plot(fit_cs_43[i].plotx .* 100, fit_cs_43[i].ploty, "C1")
    text(9, 0.42, "\$\\tau_{4, 4}=$(fit_cs_44[i].uncs[2])ms\$", color="C0")
    text(9, 0.31, "\$\\tau_{4, 3}=$(fit_cs_43[i].uncs[2] * 1000)\\mu s\$", color="C1")
    text(9, 0.18, "\$r=$(fit_cs_44[i].uncs[2] / fit_cs_43[i].uncs[2])\$")
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
    NaCsPlot.plot_survival_data(data_cs_44[i], fmt="C$(i - 1).", label="\$B_$i\$")
    plot(fit_cs_44[i].plotx, fit_cs_44[i].ploty, "C$(i - 1)")
end
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 45])
title("Cs depumping 4, 4")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_44")

figure()
for i in 1:3
    NaCsPlot.plot_survival_data(data_cs_43[i], fmt="C$(i - 1).", label="\$B_$i\$")
    plot(fit_cs_43[i].plotx, fit_cs_43[i].ploty, "C$(i - 1)")
end
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 0.55])
title("Cs depumping 4, 3")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_43")

for i in 1:3
    figure()
    NaCsPlot.plot_survival_data(data_na_22[i], fmt="C0.", label="2, 2")
    plot(fit_na_22[i].plotx, fit_na_22[i].ploty, "C0")
    NaCsPlot.plot_survival_data(NaCsData.map_params((i, v)->v * 100, data_na_21[i]),
                                fmt="C1.", label="2, 1\${}_{t\\times100}\$")
    plot(fit_na_21[i].plotx .* 100, fit_na_21[i].ploty, "C1")
    text(9, 0.42, "\$\\tau_{2, 2}=$(fit_na_22[i].uncs[2])ms\$", color="C0")
    text(9, 0.31, "\$\\tau_{2, 1}=$(fit_na_21[i].uncs[2] * 1000)\\mu s\$", color="C1")
    text(9, 0.18, "\$r=$(fit_na_22[i].uncs[2] / fit_na_21[i].uncs[2])\$")
    legend(fontsize="small")
    grid()
    ylim([0, 1])
    xlim([0, 55])
    title("Na depumping \$B_$i\$")
    xlabel("Time (ms)")
    ylabel("Survival")
    NaCsPlot.maybe_save("$(prefix)_na_b$i")
end

figure()
for i in 1:3
    NaCsPlot.plot_survival_data(data_na_22[i], fmt="C$(i - 1).", label="\$B_$i\$")
    plot(fit_na_22[i].plotx, fit_na_22[i].ploty, "C$(i - 1)")
end
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 22])
title("Na depumping 2, 2")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_22")

figure()
for i in 1:3
    NaCsPlot.plot_survival_data(data_na_21[i], fmt="C$(i - 1).", label="\$B_$i\$")
    plot(fit_na_21[i].plotx, fit_na_21[i].ploty, "C$(i - 1)")
end
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 0.55])
title("Na depumping 2, 1")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_21")

NaCsPlot.maybe_show()
