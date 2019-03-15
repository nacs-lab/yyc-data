#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190311_223739.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 ]
const specs = [(298.0563 .+ (-4:0.5:4) .* 1e-3,
                298.0563 .+ (-4:0.5:4) .* 1e-3)]
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

model_lorentzian(x, p) = p[1] .- p[2] ./ (1 .+ 4 .* ((x .- p[3]) ./ p[4]).^2)

data_05 = datas_nacs[1][1]
data_03 = datas_nacs[1][2]

fit_05 = fit_survival(model_lorentzian, data_05, [0.78, 0.2, 298.0565, 0.0015])
fit_03 = fit_survival(model_lorentzian, data_03, [0.78, 0.08, 298.056, 0.003])

const prefix = joinpath(@__DIR__, "imgs", "data_20190311_223739_short_raman")

figure()
NaCsPlot.plot_survival_data(data_05, fmt="C0.", label="0.5 ms")
plot(fit_05.plotx, fit_05.ploty, "C0-")
NaCsPlot.plot_survival_data(data_03, fmt="C1.", label="0.3 ms")
plot(fit_03.plotx, fit_03.ploty, "C1-")
text(298.0517, 0.615, ("\$f_0=$(fit_03.uncs[3])\\ \\mathrm{MHz}\$\n" *
                      "\$A=$(fit_03.uncs[2])\$\n" *
                      "\$\\Gamma=$(fit_03.uncs[4] * 1000)\\ \\mathrm{kHz}\$"),
     fontsize="small", color="C1")
text(298.0610, 0.552, ("\$A=$(fit_05.uncs[2])\$\n" *
                      "\$\\Gamma=$(fit_05.uncs[4] * 1000)\\ \\mathrm{kHz}\$\n" *
                      "\$f_0=$(fit_05.uncs[3])\\ \\mathrm{MHz}\$"),
     horizontalalignment="right",
     verticalalignment="bottom",
     fontsize="small", color="C0")
legend(fontsize="small")
grid()
ylim([0.55, 0.8])
title("Short Raman pulse")
xlabel("Frequency (MHz)")
ylabel("2-Body Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
