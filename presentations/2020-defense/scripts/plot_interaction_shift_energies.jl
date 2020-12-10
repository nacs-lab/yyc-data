#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
# using DataStructures
# using LsqFit
using LibArchive

using DelimitedFiles

function read_csv_compressed(fname, args...; kwargs...)
    LibArchive.Reader(fname) do reader
        LibArchive.support_format_raw(reader)
        LibArchive.support_filter_all(reader)
        LibArchive.next_header(reader)
        readdlm(reader, args...; kwargs...)
    end
end

const data = read_csv_compressed(joinpath(@__DIR__, "../../thesis/data/InteractionShiftLevels.csv.zst"), ',', Float64)
const idx_a0 = first(searchsorted(data[:, 1], 0))
const perturb_slope = (data[idx_a0, 2] - data[idx_a0 - 1, 2]) /
    (data[idx_a0, 1] - data[idx_a0 - 1, 1])

const prefix = joinpath(@__DIR__, "../figures/interaction_shift_energies")

function interpolate_data(a, idx)
    idx_a = first(searchsorted(data[:, 1], a))
    v1 = data[idx_a - 1, idx]
    v2 = data[idx_a, idx]
    a1 = data[idx_a - 1, 1]
    a2 = data[idx_a, 1]
    return (v1 * (a2 - a) + v2 * (a - a1)) / (a2 - a1)
end

figure(figsize=[5.6, 4.9])
plot([-0.5, 0.5], [-0.5 * perturb_slope, 0.5 * perturb_slope],
     ls="dotted", color="gray")
plot(data[:, 1], data[:, 2], "C0") # 0, 0
text(0.17, 4.7,
     "\$|m_{\\mathrm{rel},z},\\!m_{\\mathrm{COM},z}\\rangle\$\n\$=\\!|0,\\!0\\rangle\$",
     color="C0", fontsize=14)
plot(data[:, 1], data[:, 5], "C3") # 2, 0
text(-0.480, -7.1, "\$|0,\\!2\\rangle\$", color="C3", fontsize=14)
plot(data[:, 1], data[:, 6], "C4") # 1, 1
text(-0.490, 23.0, "\$|2,\\!0\\rangle\$", color="C4", fontsize=14)
plot(data[:, 1], data[:, 7], "C5") # 0, 2
text(-0.490, 39.6, "\$|1,\\!1\\rangle\$", color="C5", fontsize=14)
grid()
xlim([-0.5, 0.5])
ylim([-50, 50])
ylabel("Energies (kHz)")
xlabel("\$a / \\beta_{\\mathrm{rel},z}\$")
NaCsPlot.maybe_save("$(prefix)")

const state_prop = Dict(:ha=>"center", :va=>"center")
const ptrprops = Dict(:width=>1, :headlength=>10, :headwidth=>5, :color=>"m")

figure(figsize=[5.6, 4.9])
xinit = 0.0102
xfinal = -0.2323
yinit = interpolate_data(xinit, 2)
yfinal1 = interpolate_data(xfinal, 2)
yfinal2 = interpolate_data(xfinal, 6)
axvline(xinit, color="green", ls="dotted")
axvline(xfinal, color="red", ls="dotted")
plot(data[:, 1], data[:, 2], "C0") # 0, 0
text(0.17, 4.7,
     "\$|m_{\\mathrm{rel},z},\\!m_{\\mathrm{COM},z}\\rangle\$\n\$=\\!|0,\\!0\\rangle\$",
     color="C0", fontsize=14)
plot(data[:, 1], data[:, 5], "C3") # 2, 0
text(-0.480, -7.1, "\$|0,\\!2\\rangle\$", color="C3", fontsize=14)
plot(data[:, 1], data[:, 6], "C4") # 1, 1
text(-0.490, 23.0, "\$|2,\\!0\\rangle\$", color="C4", fontsize=14)
plot(data[:, 1], data[:, 7], "C5") # 0, 2
text(-0.490, 39.6, "\$|1,\\!1\\rangle\$", color="C5", fontsize=14)
annotate("", xy=(xfinal, yfinal1), xytext=(xinit, yinit), arrowprops=ptrprops)
annotate("", xy=(xfinal, yfinal2), xytext=(xinit, yinit), arrowprops=ptrprops)
text(xinit, -26, "Initial Spin State\nNa(2,2) Cs(4,4)",
     state_prop, color="green", rotation=90, fontsize=12, linespacing=1.5)
text(xfinal, -6, "Final Spin State\n  Na(2,2) Cs(3,3)",
     state_prop, color="red", rotation=90, fontsize=12, linespacing=1.5)
grid()
xlim([-0.5, 0.5])
ylim([-50, 50])
ylabel("Energies (kHz)")
xlabel("\$a / \\beta_{\\mathrm{rel},z}\$")
NaCsPlot.maybe_save("$(prefix)_42_32")

figure(figsize=[5.6, 4.9])
xinit = -0.2323
xfinal = 0.0046
yinit = interpolate_data(xinit, 2)
yfinal1 = interpolate_data(xfinal, 2)
yfinal2 = interpolate_data(xfinal, 6)
axvline(xinit, color="green", ls="dotted")
axvline(xfinal, color="red", ls="dotted")
plot(data[:, 1], data[:, 2], "C0") # 0, 0
text(0.17, 4.7,
     "\$|m_{\\mathrm{rel},z},\\!m_{\\mathrm{COM},z}\\rangle\$\n\$=\\!|0,\\!0\\rangle\$",
     color="C0", fontsize=14)
plot(data[:, 1], data[:, 5], "C3") # 2, 0
text(-0.480, -7.1, "\$|0,\\!2\\rangle\$", color="C3", fontsize=14)
plot(data[:, 1], data[:, 6], "C4") # 1, 1
text(-0.490, 23.0, "\$|2,\\!0\\rangle\$", color="C4", fontsize=14)
plot(data[:, 1], data[:, 7], "C5") # 0, 2
text(-0.490, 39.6, "\$|1,\\!1\\rangle\$", color="C5", fontsize=14)
annotate("", xy=(xfinal, yfinal1), xytext=(xinit, yinit), arrowprops=ptrprops)
annotate("", xy=(xfinal, yfinal2), xytext=(xinit, yinit), arrowprops=ptrprops)
text(xinit, -6, "Initial Spin State\n      Na(2,2) Cs(3,3)",
     state_prop, color="green", rotation=90, fontsize=12, linespacing=1.5)
text(xfinal, -26, "Final Spin State\nNa(1,1) Cs(3,3)",
     state_prop, color="red", rotation=90, fontsize=12, linespacing=1.5)
grid()
xlim([-0.5, 0.5])
ylim([-50, 50])
ylabel("Energies (kHz)")
xlabel("\$a / \\beta_{\\mathrm{rel},z}\$")
NaCsPlot.maybe_save("$(prefix)_32_31")

NaCsPlot.maybe_show()
