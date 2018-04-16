#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using MAT

matopen(joinpath(@__DIR__, "AT_690peak_1.mat")) do f
    global xplot, yplot, yerrplot
    xplot = read(f, "xplot")
    yplot = read(f, "yplot")
    yerrplot = read(f, "yerrplot")
end

const prefix = joinpath(@__DIR__, "eit")

figure()
errorbar(xplot[3:end], yplot[3:end], yerrplot[3:end])
xlabel("\$\\delta\$, Two photon detuning (MHz)")
ylabel("Two body survival")
ylim([0, 1])
grid()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
