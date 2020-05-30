#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

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

const fname = joinpath(@__DIR__, "20MHz_linewidth_c3SigmaOnly_3322_80kHzConfinement.csv.zst")
const data = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    readdlm(reader, ',', skipstart=1)
end
const prefix = joinpath(@__DIR__, "imgs", "raman")

figure()
l1 = plot(data[:, 1] .- 288625.081, abs.(data[:, 2] ./ 2π / 1000), "C0",
          label="\$\\Omega_{Raman}\$")
l2 = plot(data[:, 1] .- 288625.081, abs.(data[:, 4] ./ 2π / 1000), "C1",
          label="\$\\Gamma_{scatter}\$")
ylabel("\$2\\pi\\cdot \\mathrm{kHz}\$")
xlabel("One-Photon Detuning (GHz)")
ylim([0, 30])
grid()
ax = gca()
ax2 = ax.twinx()
l3 = ax2.plot(data[:, 1] .- 288625.081, abs.(data[:, 2] ./ data[:, 4]), "C2",
              label="\$\\frac{\\Omega_{Raman}}{\\Gamma_{scatter}}\$", linewidth=4)
ls = [l1; l2; l3]
legend(ls, [l.get_label() for l in ls], fontsize="small",
       loc="lower right", bbox_to_anchor=(1, 0.05))
xlim([-45, 45])
ylim([0, 60])
ax2.tick_params(axis="y", labelcolor="C2")
setp(ax2.get_yticklabels(), fontweight="bold")
NaCsPlot.maybe_save("$(prefix)_v0")

figure()
l1 = plot(data[:, 1] .- 351271.53, abs.(data[:, 2] ./ 2π / 1000), "C0",
          label="\$\\Omega_{Raman}\$")
l2 = plot(data[:, 1] .- 351271.53, abs.(data[:, 4] ./ 2π / 1000), "C1",
          label="\$\\Gamma_{scatter}\$")
ylabel("\$2\\pi\\cdot \\mathrm{kHz}\$")
xlabel("One-Photon Detuning (GHz)")
ylim([0, 800])
grid()
ax = gca()
ax2 = ax.twinx()
l3 = ax2.plot(data[:, 1] .- 351271.53, abs.(data[:, 2] ./ data[:, 4]), "C2",
              label="\$\\frac{\\Omega_{Raman}}{\\Gamma_{scatter}}\$", linewidth=4)
ls = [l1; l2; l3]
legend(ls, [l.get_label() for l in ls], fontsize="small", loc="upper right")
xlim([-12.5, 12.5])
ylim([0, 2])
ax2.tick_params(axis="y", labelcolor="C2")
setp(ax2.get_yticklabels(), fontweight="bold")
NaCsPlot.maybe_save("$(prefix)_vh")

NaCsPlot.maybe_show()
