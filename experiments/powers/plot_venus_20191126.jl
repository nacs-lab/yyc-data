#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname = joinpath(@__DIR__, "data", "venus_20191126.csv")
const data = readdlm(iname, ',', Float64, skipstart=1)

const prefix = joinpath(@__DIR__, "imgs", "venus_20191126")

figure()
plot(data[:, 1], data[:, 2], "C0o-")
grid()
ylim([0, 5])
xlabel("Wavelength (nm)")
ylabel("Output Power (W)")
tight_layout()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
