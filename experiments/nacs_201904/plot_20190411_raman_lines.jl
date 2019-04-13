#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

f = [390.0, 88, -388, -694, 600, 215, 28, -64, -115, -90, -50, 4] .+ 0.05
r = [298.191, 298.017, 297.968, 297.960, 298.95, 298.093, 298.006, 298.035, 297.999,
     298.013, 298.050, 297.985]
uncs = [0.003, 0.002, 0.002, 0.001, 0.01, 0.003, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001]

model(x, p) = p[1] .- p[2] ./ (x .- 696.5) .- p[3] ./ (x .- p[4])
fit = curve_fit(model, f, r, [297.98, 100, 10, -25])
plotx1 = linspace(-750, fit.param[4] - 5, 1000)
plotx2 = linspace(fit.param[4] + 5, 650, 1000)
ploty1 = model(plotx1, fit.param)
ploty2 = model(plotx2, fit.param)
fit_uncs = Unc.(fit.param, estimate_errors(fit), Sci)

const prefix = joinpath(@__DIR__, "imgs", "data_20190411_raman_lines")

figure()
errorbar(f, (r .- 298) .* 1000, uncs .* 1000, fmt="C0.")
plot(plotx1, (ploty1 .- 298) .* 1000, "C0")
plot(plotx2, (ploty2 .- 298) .* 1000, "C0")
text(-700, 700, ("\$f_{Raman}=f_{Raman0}-\\dfrac{a}{f_{PA}-696.5 \\mathrm{GHz}}\$\n" *
                 "   \$\\ \\ -\\dfrac{b}{f_{PA} - f_{PA0}'}\$"),
     fontsize="small")
text(-650, 300, ("\$f_{Raman0}=$(fit_uncs[1]) \\mathrm{MHz}\$\n" *
                 "\$a=$(fit_uncs[2]) \\mathrm{MHz\\cdot GHz}\$\n" *
                 "\$f_{PA0}'=$(fit_uncs[4]) \\mathrm{GHz}\$\n" *
                 "\$b=$(fit_uncs[3]) \\mathrm{MHz\\cdot GHz}\$"),
     fontsize="small")
ylim([-100, 1000])
grid()
title("Raman light shift (10 mW)")
xlabel("288XXX GHz")
ylabel("298XXX kHz")
NaCsPlot.maybe_save("$(prefix)")

figure()
errorbar(f, (r .- 298) .* 1000, uncs .* 1000, fmt="C0.")
plot(plotx1, (ploty1 .- 298) .* 1000, "C0")
plot(plotx2, (ploty2 .- 298) .* 1000, "C0")
xlim([-125, 100])
ylim([-25, 75])
grid()
title("Raman light shift (10 mW)")
xlabel("288XXX GHz")
ylabel("298XXX kHz")
NaCsPlot.maybe_save("$(prefix)_zoom")

NaCsPlot.maybe_show()
