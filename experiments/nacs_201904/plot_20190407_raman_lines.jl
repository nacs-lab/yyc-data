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
δ = r .- r0
uncs = [0.003, 0.002, 0.002, 0.001]

model(x, p) = p[1] .- p[2] ./ (x .- 698)
fit = curve_fit(model, f, δ, [0.11, 100])
plotx = linspace(-700, 450, 1000);
ploty = model(plotx, fit.param);

const prefix = joinpath(@__DIR__, "imgs", "data_20190407_raman_lines")

figure()
errorbar(f, (r .- 298) .* 1000, uncs .* 1000, fmt="C0.")
plot(plotx, (ploty .+ r0 .- 298) .* 1000, "C0")
grid()
title("Raman light shift")
xlabel("288XXX GHz")
ylabel("298XXX kHz")
NaCsPlot.maybe_save("$(prefix)_resonance")

figure()
errorbar(f, δ .* 1000, uncs .* 1000, fmt="C0.")
plot(plotx, ploty .* 1000, "C0")
grid()
title("Raman light shift")
xlabel("288XXX GHz")
ylabel("Light shift / 10mW (kHz)")
NaCsPlot.maybe_save("$(prefix)_shift")

NaCsPlot.maybe_show()
