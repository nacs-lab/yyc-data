#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190403_140224.mat",
                "data_20190403_165551.mat",
                "data_20190403_173430.mat",
                "data_20190403_184648.mat",
                "data_20190403_215913.mat",
                ]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 ]
const specs = [([0, 0.1, 0.2, 0.5, 1, 2, 5, 10],
                [0, 0.1, 0.2, 0.5, 1, 2, 5, 10]),
               ([0, 0.1, 0.2, 0.5, 1, 2, 5, 10],
                [0, 0.1, 0.2, 0.5, 1, 2, 5, 10]),
               ([0, 0.1, 0.2, 0.5, 1, 2, 5, 10] .* 20,
                [0, 0.1, 0.2, 0.5, 1, 2, 5, 10] .* 20),
               ([0.0, 20, 40, 100, 200, 500, 1000],
                [0.0, 20, 40, 100, 200, 500, 1000],
                [0.0, 10, 20, 40, 100, 200, 500, 1000],
                [500.0, 1000],
                [500.0, 1000]),
               ([0.0, 20, 40, 100, 200, 500, 1000],
                [0.0, 20, 40, 100, 200, 500, 1000],
                [0.0, 10, 20, 40, 100, 200, 500, 1000],
                [0.0, 10, 20, 40, 100, 200, 500, 1000],
                [0.0, 10, 20, 40, 100, 200, 500, 1000])]
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

const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs)

model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]

data_cs_good_lifetime = datas_cs[5][1]
data_cs_good_44_f = datas_cs[5][2]
data_cs_good_33_f = datas_cs[5][3]
data_cs_good_33_mf = datas_cs[5][4]
data_cs_good_44_mf = datas_cs[5][5]

data_cs_bad_lifetime = datas_cs[4][1]
data_cs_bad_44_f = datas_cs[4][2]
data_cs_bad_33_f = datas_cs[4][3]
data_cs_bad_33_mf = [datas_cs[2][1]; datas_cs[3][1]; datas_cs[4][4];]
data_cs_bad_44_mf = [datas_cs[2][2]; datas_cs[3][2]; datas_cs[4][5];]

const prefix = joinpath(@__DIR__, "imgs", "data_20190329_111044_flip_3322")

figure()
NaCsPlot.plot_survival_data(data_cs_bad_lifetime, fmt="C0.-", label="Bad")
NaCsPlot.plot_survival_data(data_cs_good_lifetime, fmt="C1.-", label="Good")
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 1000])
title("Lifetime")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_lifetime")

figure()
NaCsPlot.plot_survival_data(data_cs_bad_44_f, fmt="C0.-", label="Bad 4,4")
NaCsPlot.plot_survival_data(data_cs_good_44_f, fmt="C1.-", label="Good 4,4")
NaCsPlot.plot_survival_data(data_cs_bad_33_f, fmt="C2.-", label="Bad 3,3")
NaCsPlot.plot_survival_data(data_cs_good_33_f, fmt="C3.-", label="Good 3,3")
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 1000])
title("F changing")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f")

figure()
NaCsPlot.plot_survival_data(data_cs_bad_44_mf, fmt="C0.-", label="Bad 4,4")
NaCsPlot.plot_survival_data(data_cs_good_44_mf, fmt="C1.-", label="Good 4,4")
NaCsPlot.plot_survival_data(data_cs_bad_33_mf, fmt="C2.-", label="Bad 3,3")
NaCsPlot.plot_survival_data(data_cs_good_33_mf, fmt="C3.-", label="Good 3,3")
legend(fontsize="small")
grid()
ylim([0, 1])
xlim([0, 1000])
title("mF changing")
xlabel("Wait time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_mf")

NaCsPlot.maybe_show()
