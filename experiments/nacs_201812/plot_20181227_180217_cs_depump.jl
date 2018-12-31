#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit


const data_cs_align = [0.0 0.005223 0.001092
                       0.2 0.02513 0.002379
                       0.5 0.0519 0.003371
                       1.0 0.1057 0.004723
                       2.0 0.1972 0.006090
                       5.0 0.4089 0.007582
                       10 0.6345 0.007401
                       20 0.8390 0.005682]
const data_cs_mis = [0.00 0.00414 0.0009806
                     0.01 0.0779 0.004083
                     0.02 0.1401 0.005292
                     0.04 0.2590 0.006785
                     0.10 0.5134 0.007631
                     0.20 0.7401 0.006749]

function fit_survival(model, params, ratios, uncs, p0;
                      plotx=nothing, plot_lo=nothing, plot_hi=nothing, plot_scale=1.1)
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
    fit = curve_fit(model, params, ratios, uncs.^-(2/3), p0)
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

model_expm1(x, p) = p[1] .* (1 .- exp.(x ./ -p[2]))

fit_cs_align = fit_survival(model_expm1, data_cs_align[:, 1], data_cs_align[:, 2],
                            data_cs_align[:, 3], [0.9, 5])
fit_cs_mis = fit_survival(model_expm1, data_cs_mis[:, 1], data_cs_mis[:, 2],
                          data_cs_mis[:, 3], [0.9, 0.06])

const prefix = joinpath(@__DIR__, "imgs", "data_20181227_180217_cs_depump")

figure()
errorbar(data_cs_align[:, 1], data_cs_align[:, 2],
         data_cs_align[:, 3], fmt="C0.", label="\$0\\!^\\circ\$")
plot(fit_cs_align.plotx, fit_cs_align.ploty, "C0")
errorbar(data_cs_mis[:, 1] .* 100, data_cs_mis[:, 2],
         data_cs_mis[:, 3], fmt="C1.", label="\$45\\!^\\circ{}_{t\\times100}\$")
plot(fit_cs_mis.plotx .* 100, fit_cs_mis.ploty, "C1")
text(9, 0.42, "\$\\tau_{0\\!^\\circ}=$(fit_cs_align.uncs[2])\$ms", color="C0")
text(9, 0.31, "\$\\tau_{45\\!^\\circ}=$(fit_cs_mis.uncs[2])\$ms", color="C1")
text(9, 0.18, "\$r=$(fit_cs_align.uncs[2] / fit_cs_mis.uncs[2])\$")
legend()
grid()
ylim([0, 1])
xlim([0, 22])
title("Cs depumping")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
