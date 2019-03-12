#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190311_191326.mat",
                ]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 ]
const specs = [[0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 2, 4, 7, 10, 15, 25, 40, 80] .- 0.04,
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

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4)), maxcnts, specs)

model_exp(x, p) = p[1] .* exp.(x ./ -p[2])
model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]
model_exp2(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3] .* exp.(x ./ -p[4])

data_nacs = datas_nacs[1]

fit_nacs = fit_survival(model_exp2, data_nacs, [0.8, 0.4, 0.2, 70])

const prefix = joinpath(@__DIR__, "imgs", "data_20190311_191326_raman_time")

figure()
NaCsPlot.plot_survival_data(data_nacs, fmt="C0.")
plot(fit_nacs.plotx, fit_nacs.ploty, "C0")
text(0.2, 0.78, "\$\\tau_{fast}=$(fit_nacs.uncs[2])\$ ms\n\$p_{fast}=$(fit_nacs.uncs[1])\$")
text(3, 0.3, "\$\\tau_{slow}=$(fit_nacs.uncs[4])\$ ms\n\$p_{slow}=$(fit_nacs.uncs[3])\$")
grid()
ax = gca()
ax[:set_xscale]("log", nonposx="clip")
ax[:set_xticks]([0.1, 1, 10])
ax[:set_xticklabels](["0.1", "1", "10"])
ylim([0, 1])
xlim([0.05, 90])
title("Raman time")
xlabel("Time (ms)")
ylabel("2-Body Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
