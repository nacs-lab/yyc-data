#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
# using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

param_str(fit, i) = @sprintf("%.4f", fit.param[i])
params_strs(fit) =
    "\$a=$(param_str(fit, 1))\$\n\$b=$(param_str(fit, 2))\$\n\$c=$(param_str(fit, 3))\$"

const iname_f2_counterop = joinpath(@__DIR__, "data", "na_f2_counter-op_20180604.csv")
const data_f2_counterop = readdlm(iname_f2_counterop, ',', Float64, skipstart=1)
const iname_f2_coprop = joinpath(@__DIR__, "data", "na_f2_coprop_20180604.csv")
const data_f2_coprop = readdlm(iname_f2_coprop, ',', Float64, skipstart=1)
const iname_f1_coprop = joinpath(@__DIR__, "data", "na_f1_coprop_20180604.csv")
const data_f1_coprop = readdlm(iname_f1_coprop, ',', Float64, skipstart=1)
const iname_f1_up = joinpath(@__DIR__, "data", "na_f1_up_20180604.csv")
const data_f1_up = readdlm(iname_f1_up, ',', Float64, skipstart=1)
const iname_f1_down = joinpath(@__DIR__, "data", "na_f1_down_20180604.csv")
const data_f1_down = readdlm(iname_f1_down, ',', Float64, skipstart=1)

const prefix = joinpath(@__DIR__, "imgs", "na_20180604")

model(x, p) = p[1] .* sin.(p[2] .* sin.(p[3] .* x).^2)

fit_f2_counterop = curve_fit(model, data_f2_counterop[:, 1], data_f2_counterop[:, 2],
                             [maximum(data_f2_counterop[:, 2]), 2.0, 3.2])
figure()
plotmax_f2_counterop = maximum(data_f2_counterop[:, 1]) * 1.1
plotamp_f2_counterop = linspace(0, plotmax_f2_counterop, 10000)
plot(data_f2_counterop[:, 1], data_f2_counterop[:, 2], "o")
plot(plotamp_f2_counterop, model(plotamp_f2_counterop, fit_f2_counterop.param))
title("F2 Counter OP")
grid()
xlim([0, plotmax_f2_counterop])
ylim([0, ylim()[2]])
text(0.005, 2.23, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.14, 0.1, params_strs(fit_f2_counterop), fontsize=20)
xlabel("NaRaman2/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_f2_counter-op")

fit_f2_coprop = curve_fit(model, data_f2_coprop[:, 1], data_f2_coprop[:, 2],
                          [maximum(data_f2_coprop[:, 2]), 2.0, 4.3])
figure()
plotmax_f2_coprop = maximum(data_f2_coprop[:, 1]) * 1.1
plotamp_f2_coprop = linspace(0, plotmax_f2_coprop, 10000)
plot(data_f2_coprop[:, 1], data_f2_coprop[:, 2], "o")
plot(plotamp_f2_coprop, model(plotamp_f2_coprop, fit_f2_coprop.param))
title("F2 Co-prop")
grid()
xlim([0, plotmax_f2_coprop])
ylim([0, ylim()[2]])
text(0.01, 1.55, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.14, 0.1, params_strs(fit_f2_coprop), fontsize=20)
xlabel("NaRaman2/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_f2_coprop")

fit_f1_coprop = curve_fit(model, data_f1_coprop[:, 1], data_f1_coprop[:, 2],
                          [maximum(data_f1_coprop[:, 2]), 2.0, 3.9])
figure()
plotmax_f1_coprop = maximum(data_f1_coprop[:, 1]) * 1.1
plotamp_f1_coprop = linspace(0, plotmax_f1_coprop, 10000)
plot(data_f1_coprop[:, 1], data_f1_coprop[:, 2], "o")
plot(plotamp_f1_coprop, model(plotamp_f1_coprop, fit_f1_coprop.param))
title("F1 Co-prop")
grid()
xlim([0, plotmax_f1_coprop])
ylim([0, ylim()[2]])
text(0.01, 4.3, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.14, 0.1, params_strs(fit_f1_coprop), fontsize=20)
xlabel("NaRaman1/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_f1_coprop")

fit_f1_up = curve_fit(model, data_f1_up[:, 1], data_f1_up[:, 2],
                      [maximum(data_f1_up[:, 2]), 2.0, 1.6])
figure()
plotmax_f1_up = maximum(data_f1_up[:, 1]) * 1.1
plotamp_f1_up = linspace(0, plotmax_f1_up, 10000)
plot(data_f1_up[:, 1], data_f1_up[:, 2], "o")
plot(plotamp_f1_up, model(plotamp_f1_up, fit_f1_up.param))
title("F1 Up")
grid()
xlim([0, plotmax_f1_up])
ylim([0, ylim()[2]])
text(0.35, 3.0, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.41, 0.1, params_strs(fit_f1_up), fontsize=20)
xlabel("NaRaman1/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_f1_up")

fit_f1_down = curve_fit(model, data_f1_down[:, 1], data_f1_down[:, 2],
                        [maximum(data_f1_down[:, 2]), 2.0, 4.0])
figure()
plotmax_f1_down = maximum(data_f1_down[:, 1]) * 1.1
plotamp_f1_down = linspace(0, plotmax_f1_down, 10000)
plot(data_f1_down[:, 1], data_f1_down[:, 2], "o")
plot(plotamp_f1_down, model(plotamp_f1_down, fit_f1_down.param))
title("F1 Down")
grid()
xlim([0, plotmax_f1_down])
ylim([0, ylim()[2]])
text(0.01, 6.3, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.18, 0.1, params_strs(fit_f1_down), fontsize=20)
xlabel("NaRaman1/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_f1_down")

NaCsPlot.maybe_show()
