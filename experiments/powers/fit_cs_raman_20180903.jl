#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using LsqFit
import NaCsCalc.Format: Unc, Sci

param_str(fit, i) = @sprintf("%.4f", fit.param[i])
params_strs(fit) =
    "\$a=$(param_str(fit, 1))\$\n\$b=$(param_str(fit, 2))\$\n\$c=$(param_str(fit, 3))\$"

const iname_f3_counterop = joinpath(@__DIR__, "data", "cs_f3_counterop_20180903.csv")
const data_f3_counterop = readdlm(iname_f3_counterop, ',', Float64, skipstart=1)
const iname_f4_coprop = joinpath(@__DIR__, "data", "cs_f4_coprop_20180903.csv")
const data_f4_coprop = readdlm(iname_f4_coprop, ',', Float64, skipstart=1)
const iname_f3_coprop = joinpath(@__DIR__, "data", "cs_f3_coprop_20180903.csv")
const data_f3_coprop = readdlm(iname_f3_coprop, ',', Float64, skipstart=1)
const iname_f4_up = joinpath(@__DIR__, "data", "cs_f4_up_20180903.csv")
const data_f4_up = readdlm(iname_f4_up, ',', Float64, skipstart=1)
const iname_f4_down = joinpath(@__DIR__, "data", "cs_f4_down_20180903.csv")
const data_f4_down = readdlm(iname_f4_down, ',', Float64, skipstart=1)

const prefix = joinpath(@__DIR__, "imgs", "cs_20180903")

model(x, p) = p[1] .* sin.(p[2] .* sin.(p[3] .* x).^2)
# model′(x, p) = p[1] .* sin.(p[2] .* x.^2)

fit_f3_counterop = curve_fit(model, data_f3_counterop[:, 1], data_f3_counterop[:, 2],
                             [maximum(data_f3_counterop[:, 2]), 1.5, 1.12])
figure()
plotmax_f3_counterop = maximum(data_f3_counterop[:, 1]) * 1.1
plotamp_f3_counterop = linspace(0, plotmax_f3_counterop, 10000)
plot(data_f3_counterop[:, 1], data_f3_counterop[:, 2], "o")
plot(plotamp_f3_counterop, model(plotamp_f3_counterop, fit_f3_counterop.param))
title("F3 Counter OP")
grid()
xlim([0, plotmax_f3_counterop])
ylim([0, ylim()[2]])
text(0.06, 15, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.65, 1, params_strs(fit_f3_counterop), fontsize=20)
xlabel("CsRamanF3/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_f3_counterop")

fit_f4_coprop = curve_fit(model, data_f4_coprop[:, 1], data_f4_coprop[:, 2],
                          [maximum(data_f4_coprop[:, 2]), 1.6, 1.0])
figure()
plotmax_f4_coprop = maximum(data_f4_coprop[:, 1]) * 1.1
plotamp_f4_coprop = linspace(0, plotmax_f4_coprop, 10000)
plot(data_f4_coprop[:, 1], data_f4_coprop[:, 2], "o")
plot(plotamp_f4_coprop, model(plotamp_f4_coprop, fit_f4_coprop.param))
title("F4 Co-prop")
grid()
xlim([0, plotmax_f4_coprop])
ylim([0, ylim()[2]])
text(0.06, 20, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.55, 1, params_strs(fit_f4_coprop), fontsize=20)
xlabel("CsRamanF4/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_f4_coprop")

fit_f3_coprop = curve_fit(model, data_f3_coprop[:, 1], data_f3_coprop[:, 2],
                          [maximum(data_f3_coprop[:, 2]), 1.0, 1.0])
figure()
plotmax_f3_coprop = maximum(data_f3_coprop[:, 1]) * 1.1
plotamp_f3_coprop = linspace(0, plotmax_f3_coprop, 10000)
plot(data_f3_coprop[:, 1], data_f3_coprop[:, 2], "o")
plot(plotamp_f3_coprop, model(plotamp_f3_coprop, fit_f3_coprop.param))
title("F3 Co-prop")
grid()
xlim([0, plotmax_f3_coprop])
ylim([0, ylim()[2]])
text(0.06, 15, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.6, 1, params_strs(fit_f3_coprop), fontsize=20)
xlabel("CsRamanF3/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_f3_coprop")

fit_f4_up = curve_fit(model, data_f4_up[:, 1], data_f4_up[:, 2],
                      [maximum(data_f4_up[:, 2]), 1.3, 1.0])
figure()
plotmax_f4_up = maximum(data_f4_up[:, 1]) * 1.1
plotamp_f4_up = linspace(0, plotmax_f4_up, 10000)
plot(data_f4_up[:, 1], data_f4_up[:, 2], "o")
plot(plotamp_f4_up, model(plotamp_f4_up, fit_f4_up.param))
title("F4 Up")
grid()
xlim([0, plotmax_f4_up])
ylim([0, ylim()[2]])
text(0.06, 17, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.6, 1, params_strs(fit_f4_up), fontsize=20)
xlabel("CsRamanF4/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_f4_up")

fit_f4_down = curve_fit(model, data_f4_down[:, 1], data_f4_down[:, 2],
                        [maximum(data_f4_down[:, 2]), 1.0, 1.0])
figure()
plotmax_f4_down = maximum(data_f4_down[:, 1]) * 1.1
plotamp_f4_down = linspace(0, plotmax_f4_down, 10000)
plot(data_f4_down[:, 1], data_f4_down[:, 2], "o")
plot(plotamp_f4_down, model(plotamp_f4_down, fit_f4_down.param))
title("F4 Down")
grid()
xlim([0, plotmax_f4_down])
ylim([0, ylim()[2]])
text(0.05, 13, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.6, 1, params_strs(fit_f4_down), fontsize=20)
xlabel("CsRamanF4/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_f4_down")

NaCsPlot.maybe_show()
