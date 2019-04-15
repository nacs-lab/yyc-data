#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190415_103612.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 ]
const specs = [([0, 0.2, 0.5, 1, 2, 5, 10, 20, 30, 40],
                [0, 1, 2, 5, 10, 20, 30, 40] .* 10)]
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

model_exp1(x, p) = p[1] .* exp.(x ./ -p[2])
model_exp_off(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3]
model_exp2(x, p) = p[1] .* exp.(x ./ -p[2]) .+ p[3] .* exp.(x ./ -p[4])

data_nacs_cool = datas_nacs[1][1]
data_nacs_nocool = datas_nacs[1][2]

fit_nacs_cool = fit_survival(model_exp_off, data_nacs_cool, [0.75, 10.0, 0.1])
fit_nacs_nocool = fit_survival(model_exp_off, data_nacs_nocool, [0.75, 130, 0.1])

# @show fit_nacs_cool.uncs
# @show fit_nacs_nocool.uncs

const prefix = joinpath(@__DIR__, "imgs", "data_20190415_103612_pa_times")

figure()
NaCsPlot.plot_survival_data(data_nacs_cool, fmt="C0.", label="w/ RSC")
plot(fit_nacs_cool.plotx, fit_nacs_cool.ploty, "C0")
NaCsPlot.plot_survival_data(data_nacs_nocool, fmt="C1.", label="w/o RSC")
plot(fit_nacs_nocool.plotx, fit_nacs_nocool.ploty, "C1")
text(0.2, 0.5, "\$\\tau=$(fit_nacs_cool.uncs[2])\\ \\mathrm{ms}\$",
     color="C0", fontsize="small")
text(8, 0.78, "\$\\tau=$(fit_nacs_nocool.uncs[2] / 1000)\\ \\mathrm{s}\$",
     color="C1", fontsize="small")
legend(fontsize="small", loc="lower left")
grid()
ylim([0.2, 0.9])
xlim([0.1, 500])
ax = gca()
ax[:set_xscale]("log", nonposx="clip")
ax[:set_xticks]([0.1, 1, 10, 100])
ax[:set_xticklabels](["0.1", "1", "10", "100"])
title("PA rates")
xlabel("Wait time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
