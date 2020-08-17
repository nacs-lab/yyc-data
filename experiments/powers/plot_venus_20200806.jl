#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using DelimitedFiles
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname0 = joinpath(@__DIR__, "data", "venus_20191126.csv")
const data0 = readdlm(iname0, ',', Float64, skipstart=1)

const iname = joinpath(@__DIR__, "data", "venus_20200806.csv")
const data = readdlm(iname, ',', Float64, skipstart=1)

const prefix = joinpath(@__DIR__, "imgs", "venus_20200806")

figure()
plot(data[:, 1], data[:, 2], "C0o-")
grid()
ylim([0, 3])
xlabel("Wavelength (nm)")
ylabel("Output Power (W)")
tight_layout()
NaCsPlot.maybe_save("$(prefix)")

figure()
plot(data0[:, 1], data0[:, 2], "C1o-", label="2019/11/26")
plot(data[:, 1], data[:, 2], "C0o-", label="2020/08/06")
grid()
ylim([0, 5])
xlabel("Wavelength (nm)")
ylabel("Output Power (W)")
legend(fontsize="x-small")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_cmp")

NaCsPlot.maybe_show()
