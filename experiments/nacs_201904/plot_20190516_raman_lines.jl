#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

f = [447.5, 449, 452, 455, 458]
r = [298.185, 298.21, 298.285, 298.39, 298.58]
uncs = [0.005, 0.005, 0.005, 0.005, 0.007]

model(x, p) = p[1] .- p[2] ./ (x .- p[3])
fit = curve_fit(model, f, r, [298.0, 100, 470])
plotx = linspace(445, 460, 1000)
ploty = model(plotx, fit.param)
fit_uncs = Unc.(fit.param, estimate_errors(fit), Sci)

const prefix = joinpath(@__DIR__, "imgs", "data_20190516_raman_lines")

figure()
errorbar(f, (r .- 298) .* 1000, uncs .* 1000, fmt="C0.")
plot(plotx, (ploty .- 298) .* 1000, "C0")
text(445, 700, "\$f_{Raman}=f_{Raman0}-\\dfrac{a}{f_{PA}-f_{PA0}}\$",
     fontsize="small")
text(447, 400, ("\$f_{Raman0}=$(fit_uncs[1]) \\mathrm{MHz}\$\n" *
                 "\$a=$(fit_uncs[2]) \\mathrm{MHz\\cdot GHz}\$\n" *
                 "\$f_{PA0}=$(fit_uncs[3]) \\mathrm{GHz}\$"),
     fontsize="small")
grid()
title("Raman light shift (11 mW)")
xlabel("307XXX GHz")
ylabel("298XXX kHz")
NaCsPlot.maybe_save("$(prefix)_resonance")

NaCsPlot.maybe_show()
