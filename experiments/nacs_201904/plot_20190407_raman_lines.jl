#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

f = [390.0, 88, -388, -694]
r = [298.191, 298.017, 297.968, 297.960]
r0 = 297.864 * 2 - 297.960
sr0 = 0.001 * √(5)
δ = r .- r0
uncs = [0.003, 0.002, 0.002, 0.001]

model(x, p) = p[1] .- p[2] ./ (x .- 696.5)
fit = curve_fit(model, f, r, [297.98, 100])
plotx = linspace(-700, 450, 1000)
ploty = model(plotx, fit.param)
fit_uncs = Unc.(fit.param, estimate_errors(fit), Sci)

const prefix = joinpath(@__DIR__, "imgs", "data_20190407_raman_lines")

figure()
errorbar(f, (r .- 298) .* 1000, uncs .* 1000, fmt="C0.")
plot(plotx, (ploty .- 298) .* 1000, "C0")
text(-700, 200, "\$f_{Raman}=f_{Raman0}-\\dfrac{a}{f_{PA}-696.5 \\mathrm{GHz}}\$",
     fontsize="small")
text(-650, 100, ("\$f_{Raman0}=$(fit_uncs[1]) \\mathrm{MHz}\$\n" *
                 "\$a=$(fit_uncs[2]) \\mathrm{MHz\\cdot GHz}\$"),
     fontsize="small")
grid()
title("Raman light shift (10 mW)")
xlabel("288XXX GHz")
ylabel("298XXX kHz")
NaCsPlot.maybe_save("$(prefix)_resonance")

figure()
errorbar(f, δ .* 1000, uncs .* 1000, fmt="C0.")
plot(plotx, (ploty .- r0) .* 1000, "C0")
text(-700, 400, "\$f_{Raman}(P=0)=$(Unc(r0, sr0, Sci)) \\mathrm{MHz}\$",
     fontsize="small")
grid()
title("Raman light shift (10 mW)")
xlabel("288XXX GHz")
ylabel("Light shift / 10mW (kHz)")
NaCsPlot.maybe_save("$(prefix)_shift")

NaCsPlot.maybe_show()
