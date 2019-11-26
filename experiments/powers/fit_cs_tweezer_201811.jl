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
params_strs2(fit) =
    "\$a=$(param_str(fit, 1))\$\n\$b=$(param_str(fit, 2))\$"

const iname_paaom = joinpath(@__DIR__, "data", "cs_tweezer_paaom_201811.csv")
const data_paaom = readdlm(iname_paaom, ',', Float64, skipstart=1)
const iname_padpaom = joinpath(@__DIR__, "data", "cs_tweezer_padpaom_201811.csv")
const data_padpaom = readdlm(iname_padpaom, ',', Float64, skipstart=1)
const iname_padpaom2 = joinpath(@__DIR__, "data", "cs_tweezer_padpaom2_201811.csv")
const data_padpaom2 = readdlm(iname_padpaom2, ',', Float64, skipstart=1)

@assert data_padpaom2[:, 1] == data_padpaom[:, 1]
data_padpaom2[:, 2] .-= data_padpaom[:, 2]

const prefix = joinpath(@__DIR__, "imgs", "cs_tweezer_201811")

model(x, p) = p[1] .* sin.(p[2] .* sin.(p[3] .* x).^2)
model2(x, p) = p[1] .* sin.(p[2] .* x).^4
model3(x, p) = p[1] .* (1 .- sin.(p[2] .* x).^2).^2

fit_paaom = curve_fit(model, data_paaom[:, 1], data_paaom[:, 2],
                      [maximum(data_paaom[:, 2]), 1.5, 2.1])
figure()
plotmax_paaom = maximum(data_paaom[:, 1]) * 1.1
plotamp_paaom = linspace(0, plotmax_paaom, 10000)
plot(data_paaom[:, 1], data_paaom[:, 2], "o")
plot(plotamp_paaom, model(plotamp_paaom, fit_paaom.param))
title("PAAOM")
grid()
xlim([0, plotmax_paaom])
ylim([0, ylim()[2]])
text(0.35, 19, "\$a\\cdot\\sin(b\\cdot\\sin^2(c\\cdot AMP))\$", fontsize=18)
text(0.45, 3, params_strs(fit_paaom), fontsize=20)
xlabel("PAAOM/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_paaom")

fit_padpaom = curve_fit(model2, data_padpaom[:, 1], data_padpaom[:, 2],
                        [maximum(data_padpaom[:, 2]), 1.6])
figure()
plotmax_padpaom = maximum(data_padpaom[:, 1]) * 1.1
plotamp_padpaom = linspace(0, plotmax_padpaom, 10000)
plot(data_padpaom[:, 1], data_padpaom[:, 2], "o")
plot(plotamp_padpaom, model2(plotamp_padpaom, fit_padpaom.param))
title("PADPAOM 1st order")
grid()
xlim([0, plotmax_padpaom])
ylim([0, ylim()[2]])
text(0.05, 40, "\$a\\cdot\\sin^4(b\\cdot AMP)\$", fontsize=18)
text(0.08, 28, params_strs2(fit_padpaom), fontsize=20)
xlabel("PADPAOM/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_padpaom")

fit_padpaom2 = curve_fit(model3, data_padpaom2[:, 1], data_padpaom2[:, 2],
                         [maximum(data_padpaom2[:, 2]), 1.0])
figure()
plotmax_padpaom2 = maximum(data_padpaom2[:, 1]) * 1.1
plotamp_padpaom2 = linspace(0, plotmax_padpaom2, 10000)
plot(data_padpaom2[:, 1], data_padpaom2[:, 2], "o")
plot(plotamp_padpaom2, model3(plotamp_padpaom2, fit_padpaom2.param))
title("PADPAOM 0th order")
grid()
xlim([0, plotmax_padpaom2])
ylim([0, ylim()[2]])
text(0.2, 38, "\$a\\cdot(1-\\sin^2(b\\cdot AMP))^2\$", fontsize=18)
text(0.25, 28, params_strs2(fit_padpaom2), fontsize=20)
xlabel("PADPAOM/AMP")
ylabel("Power (mW)")
NaCsPlot.maybe_save("$(prefix)_padpaom2")

NaCsPlot.maybe_show()
