#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

f = [306540.7, 306560.7, 306580.7, 306530.7, 306460.7, 306480.7, 306490.7] .- 0.02
r = [367.405, 367.579, 367.6484, 367.196, 368.159, 368.439, 368.885]
uncs = [0.02, 0.045, 0.0066, 0.022, 0.01, 0.012, 0.013]

model(x, p) = p[1] .- p[2] ./ (x .- p[3])
fit = curve_fit(model, f, r, [367.8, 12, 306510.0])
plotx1 = linspace(306450, 306496, 1000)
ploty1 = model(plotx1, fit.param)
plotx2 = linspace(306520, 306600, 1000)
ploty2 = model(plotx2, fit.param)
fit_uncs = Unc.(fit.param, estimate_errors(fit), Sci)

# @show fit_uncs

const prefix = joinpath(@__DIR__, "imgs", "data_20191226_raman_lines")

figure()
errorbar(f .- 306000, r, uncs, fmt="C0.")
plot(plotx1 .- 306000, ploty1, "C0")
plot(plotx2 .- 306000, ploty2, "C0")
text(450, 367.72, "\$f_{Raman}=f_{Raman0}-\\dfrac{a}{f_{PA}-f_{PA0}}\$",
     fontsize="small", color="C0")
text(500, 368.66, ("\$f_{Raman0}=$(fit_uncs[1]) \\mathrm{MHz}\$\n" *
                   "\$a=$(fit_uncs[2]) \\mathrm{MHz\\cdot GHz}\$\n" *
                   "\$f_{PA0}=$(fit_uncs[3]) \\mathrm{GHz}\$"),
     fontsize="small", color="C0")
grid()
title("Raman light shift (8 mW)")
# ylim([-1100, 800])
xlabel("Up leg frequency (306XXX GHz)")
ylabel("Raman resonance (MHz)")
NaCsPlot.maybe_save("$(prefix)_resonance")

NaCsPlot.maybe_show()
