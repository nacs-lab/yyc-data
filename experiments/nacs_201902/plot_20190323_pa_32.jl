#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190323_185612.mat",
                "data_20190324_093902.mat",
                "data_20190324_153235.mat",
                ]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 ]
const specs = [[0, 0.2, 0.5, 1, 2, 4, 6, 8, 10],
               [0, 0.2, 0.5, 1, 2, 4, 6, 8, 10] .* 0.5,
               [0, 0.2, 0.5, 1, 2, 4, 6, 8, 10] .* 0.5,
               ]
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

const datas_cs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]

data_620 = datas_cs[1]
data_680 = datas_cs[2]
data_690 = datas_cs[3]

fit_620 = fit_survival(model_exp_off, data_620, [0.5, 2, 0.5])
fit_680 = fit_survival(model_exp_off, data_680, [0.5, 0.5, 0.5])
fit_690 = fit_survival(model_exp_off, data_690, [0.5, 0.5, 0.5])

const prefix = joinpath(@__DIR__, "imgs", "data_20190323_pa_32")

figure()
NaCsPlot.plot_survival_data(data_620, fmt="C0.")
plot(fit_620.plotx, fit_620.ploty, "C0-", label="620 GHz")
NaCsPlot.plot_survival_data(data_680, fmt="C1.")
plot(fit_680.plotx, fit_680.ploty, "C1-", label="680 GHz")
NaCsPlot.plot_survival_data(data_690, fmt="C2.")
plot(fit_690.plotx, fit_690.ploty, "C2-", label="690 GHz")
text(2, 0.47, "\$\\tau_{620}=$(fit_620.uncs[2]) ms\$", color="C0")
text(3.5, 0.36, "\$\\tau_{680}=$(fit_680.uncs[2]) ms\$", color="C1")
text(0.8, 0.31, "\$\\tau_{690}=$(fit_690.uncs[2]) ms\$", color="C2")
legend(fontsize="small")
grid()
# ax = gca()
# ax[:set_xscale]("log", nonposx="clip")
# ax[:set_xticks]([0.1, 1, 10])
# ax[:set_xticklabels](["0.1", "1", "10"])
ylim([0.30, 0.65])
xlim([0, 7])
title("3, 3 + 2, 2 loss")
xlabel("Wait time (ms)")
ylabel("Two body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
