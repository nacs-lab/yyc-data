#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

f = [447.5, 449, 452, 455, 458, 475]
r = [298.185, 298.21, 298.285, 298.39, 298.58, 297.07]
uncs = [0.005, 0.005, 0.005, 0.005, 0.007, 0.005]

model(x, p) = p[1] .- p[2] ./ (x .- p[3])
fit = curve_fit(model, f, r, [298.0, 100, 466])
plotx1 = linspace(445, 464, 1000)
ploty1 = model(plotx1, fit.param)
plotx2 = linspace(470, 480, 1000)
ploty2 = model(plotx2, fit.param)
fit_uncs = Unc.(fit.param, estimate_errors(fit), Sci)

const prefix = joinpath(@__DIR__, "imgs", "data_20190516_raman_lines")

figure()
errorbar(f, (r .- 298) .* 1000, uncs .* 1000, fmt="C0.")
plot(plotx1, (ploty1 .- 298) .* 1000, "C0")
plot(plotx2, (ploty2 .- 298) .* 1000, "C0")
text(450, -100, "\$f_{Raman}=f_{Raman0}-\\dfrac{a}{f_{PA}-f_{PA0}}\$",
     fontsize="small")
text(452, -700, ("\$f_{Raman0}=$(fit_uncs[1]) \\mathrm{MHz}\$\n" *
                 "\$a=$(fit_uncs[2]) \\mathrm{MHz\\cdot GHz}\$\n" *
                 "\$f_{PA0}=$(fit_uncs[3]) \\mathrm{GHz}\$"),
     fontsize="small")
grid()
title("Raman light shift (11 mW)")
xlim([445, 480])
ylim([-1100, 800])
xlabel("307XXX GHz")
ylabel("298XXX kHz")
NaCsPlot.maybe_save("$(prefix)_resonance")

NaCsPlot.maybe_show()
