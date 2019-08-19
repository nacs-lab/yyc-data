#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using LsqFit

const powers = [3, 7.5, 15]
const linewidths = [(76 / 45^2 + 70 / 19^2) / (1 / 45^2 + 1 / 19^2),
                    152,
                    (263 / 76^2 + 380 / 130^2) / (1 / 76^2 + 1 / 130^2)]
const uncs = [1 / sqrt(1 / 45^2 + 1 / 19^2),
              26,
              1 / sqrt(1 / 76^2 + 1 / 130^2)]

const prefix = joinpath(@__DIR__, "imgs", "data_20190819_pa_linewidths")

model_lin0(x, p) = p[1] .* x
model_sq(x, p) = sqrt.(p[1]^2 .+ (p[2] .* x).^4)

fit_lin0 = curve_fit(model_lin0, powers, linewidths, uncs.^-(2/3), [20.0])
fit_sq = curve_fit(model_sq, powers, linewidths, uncs.^-(2/3), [50, 2.0])

fit_unc_lin0 = Unc.(fit_lin0.param, estimate_errors(fit_lin0), Sci)

plotx = linspace(0, 16, 1000)

figure()
errorbar(powers, linewidths, uncs, fmt="C0.")
plot(plotx, model_lin0(plotx, fit_lin0.param), "C0")
plot(plotx, model_sq(plotx, fit_sq.param), "C1")
axhline(50, color="C2", ls="--")
text(0.6, 220, "$(fit_unc_lin0[1]) MHz/mW", color="C0")
xlabel("Tweezer Power (mW)")
ylabel("Linewidth (MHz)")
title("PA linewidth")
xlim([0, 16])
ylim([0, 380])
grid()
tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
