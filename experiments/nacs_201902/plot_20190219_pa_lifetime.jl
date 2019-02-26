#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190219_151703.mat",
                "data_20190220_094017.mat",
                "data_20190220_144236.mat",
                "data_20190220_153441.mat",
                "data_20190220_171238.mat",
                "data_20190220_215125.mat", # 682 GHz
                "data_20190221_081927.mat",
                "data_20190221_143511.mat",
                "data_20190221_154307.mat",
                "data_20190221_223037.mat",
                "data_20190222_110328.mat",
                "data_20190222_185937.mat"
                ]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int), # 682 GHz
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int)
                 ]
const specs = [([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500], # 682 GHz
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 10, 20, 50, 100, 200, 500],
                [0, 10, 20, 50, 100, 200, 500]),
               ([0, 5, 10, 20, 50, 100, 200],
                [0, 5, 10, 20, 50, 100, 200]),
               [0, 1, 2, 4, 6, 8, 10, 12, 15, 20]]
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

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4)), maxcnts, specs)

model_exp(x, p) = p[1] .* exp.(x ./ -p[2])
model_exp2(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3] .* exp.(x ./ -p[4])

data_640 = datas_nacs[1]
data_656 = datas_nacs[2]
data_670 = [[datas_nacs[3][1]; datas_nacs[4][1]; datas_nacs[5][1]],
            [datas_nacs[3][2]; datas_nacs[4][2]; datas_nacs[5][2]]]
data_682 = datas_nacs[6]
data_688 = [[datas_nacs[7][1]; datas_nacs[8][1]; datas_nacs[9][1]],
            [datas_nacs[7][2]; datas_nacs[8][2]; datas_nacs[9][2]]]
data_690 = datas_nacs[10]
data_691_5 = [datas_nacs[11][1], [datas_nacs[11][2]; datas_nacs[12]]]

fit_640_hot = fit_survival(model_exp, data_640[2], [0.75, 800])
fit_640_cold = fit_survival(model_exp2, data_640[4], [0.5, 200, 0.2, 800])
fit_656_hot = fit_survival(model_exp, data_656[1], [0.75, 800])
fit_656_cold = fit_survival(model_exp2, data_656[2], [0.5, 200, 0.2, 800])
fit_670_hot = fit_survival(model_exp, data_670[1], [0.75, 800])
fit_670_cold = fit_survival(model_exp2, data_670[2], [0.5, 200, 0.2, 800])
fit_682_hot = fit_survival(model_exp, data_682[1], [0.75, 800])
fit_682_cold = fit_survival(model_exp2, data_682[2], [0.5, 200, 0.2, 800])
fit_688_hot = fit_survival(model_exp, data_688[1], [0.75, 800])
fit_688_cold = fit_survival(model_exp2, data_688[2], [0.5, 200, 0.2, 800])
fit_690_hot = fit_survival(model_exp, data_690[1], [0.75, 800])
fit_690_cold = fit_survival(model_exp2, data_690[2], [0.5, 200, 0.2, 800])
fit_691_5_hot = fit_survival(model_exp, data_691_5[1], [0.75, 800])
fit_691_5_cold = fit_survival(model_exp2, data_691_5[2], [0.5, 200, 0.2, 800])

const prefix = joinpath(@__DIR__, "imgs", "data_20190219_ps_lifetime")

figure()
NaCsPlot.plot_survival_data(data_640[2], fmt="C0.", label="640")
plot(fit_640_hot.plotx, fit_640_hot.ploty, "C0")
NaCsPlot.plot_survival_data(data_656[1], fmt="C1.", label="656")
plot(fit_656_hot.plotx, fit_656_hot.ploty, "C1")
NaCsPlot.plot_survival_data(data_670[1], fmt="C2.", label="670")
plot(fit_670_hot.plotx, fit_670_hot.ploty, "C2")
NaCsPlot.plot_survival_data(data_682[1], fmt="C3.", label="682")
plot(fit_682_hot.plotx, fit_682_hot.ploty, "C3")
NaCsPlot.plot_survival_data(data_688[1], fmt="C4.", label="688")
plot(fit_688_hot.plotx, fit_688_hot.ploty, "C4")
NaCsPlot.plot_survival_data(data_690[1], fmt="C5.", label="690")
plot(fit_690_hot.plotx, fit_690_hot.ploty, "C5")
NaCsPlot.plot_survival_data(data_691_5[1], fmt="C6.", label="691.5")
plot(fit_691_5_hot.plotx, fit_691_5_hot.ploty, "C6")
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.2)
grid()
ylim([0, 1])
xlim([0, 550])
title("2-Body Lifetime (Hot)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_hot1")

figure()
NaCsPlot.plot_survival_data(data_640[4], fmt="C0.", label="640")
plot(fit_640_cold.plotx, fit_640_cold.ploty, "C0")
NaCsPlot.plot_survival_data(data_656[2], fmt="C1.", label="656")
plot(fit_656_cold.plotx, fit_656_cold.ploty, "C1")
NaCsPlot.plot_survival_data(data_670[2], fmt="C2.", label="670")
plot(fit_670_cold.plotx, fit_670_cold.ploty, "C2")
NaCsPlot.plot_survival_data(data_682[2], fmt="C3.", label="682")
plot(fit_682_cold.plotx, fit_682_cold.ploty, "C3")
NaCsPlot.plot_survival_data(data_688[2], fmt="C4.", label="688")
plot(fit_688_cold.plotx, fit_688_cold.ploty, "C4")
NaCsPlot.plot_survival_data(data_690[2], fmt="C5.", label="690")
plot(fit_690_cold.plotx, fit_690_cold.ploty, "C5")
NaCsPlot.plot_survival_data(data_691_5[2], fmt="C6.", label="691.5")
plot(fit_691_5_cold.plotx, fit_691_5_cold.ploty, "C6")
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.2)
grid()
ylim([0, 1])
xlim([0, 250])
title("2-Body Lifetime (Cold)")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cold1")

NaCsPlot.maybe_show()
