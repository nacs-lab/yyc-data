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

read_zst_csv(fname) = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    readdlm(reader, ',', skipstart=1)
end

const fname = joinpath(@__DIR__, "../../molecular_raman_paper/data/50MHz_linewidth_8TotalExcitedStates_3322_80kHzConfinement_3.75mWPowerAtAtom.csv.zst")
const fname2 = joinpath(@__DIR__, "../../molecular_raman_paper/data/50MHz_linewidth_8TotalExcitedStates_3322_80kHzConfinement_3.75mWPowerAtAtom_ThresholdContribution.csv.zst")
const data = read_zst_csv(fname)
const data2 = read_zst_csv(fname2)
const prefix = joinpath(@__DIR__, "../figures/raman_transfer_theory_ratio")

figure()
ax1 = gca()
plot(data[:, 1] .- 339724.57, abs.(data[:, 2] ./ 2π / 1000), "C1",
     label="\$v'=40\$")
plot(data[:, 1] .- 288625.081, abs.(data[:, 2] ./ 2π / 1000), "C0",
     label="\$v'=0\$")
grid()
yscale("log")
ax1.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
xlim([-300, 300])
ylim([6, 20000])
xlabel("One-Photon Detuning (GHz)")
ylabel("\$\\Omega_{\\mathrm{Raman}}~(2\\pi\\!\\times\\!\\mathrm{kHz})\$")
legend(fontsize=13.88, loc="center left", handlelength=1, handletextpad=0.3)
NaCsPlot.maybe_save("$(prefix)_omega")

# figure()
# ax2 = gca()
# plot(data[:, 1] .- 339724.57, abs.(data[:, 4] ./ 2π / 1000), "C1",
#      label="\$v'=40\$")
# plot(data[:, 1] .- 288625.081, abs.(data[:, 4] ./ 2π / 1000), "C0",
#      label="\$v'=0\$")
# grid()
# yscale("log")
# ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# xlim([-300, 300])
# ylim([0.12, 200000])
# xlabel("One-Photon Detuning (GHz)")
# ylabel("\$\\Gamma_{\\mathrm{scatter}}~(2\\pi\\!\\times\\!\\mathrm{kHz})\$")
# legend(fontsize=13.88, loc="upper left", bbox_to_anchor=(0, 0.8),
#        handlelength=1, handletextpad=0.3)
# NaCsPlot.maybe_save("$(prefix)_gamma")

figure()
ax3 = gca()
plot(data[:, 1] .- 288625.081, abs.(data[:, 2] ./ data[:, 4]), "C0",
     label="\$v'=0\$")
plot(data[:, 1] .- 339724.57, abs.(data[:, 2] ./ data[:, 4]), "C1",
     label="\$v'=40\$")
grid()
yscale("log")
ax3.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
yticks([0.1, 1, 10], ["0.1", "1", "10"])
xlim([-300, 300])
ylim([0.02, 30])
xlabel("One-Photon Detuning (GHz)")
ylabel("\$\\Omega_{\\mathrm{Raman}}/\\Gamma_{\\mathrm{scatter}}\$")
legend(fontsize=13.88, loc="upper left", bbox_to_anchor=(0, 0.8),
       handlelength=1, handletextpad=0.3)
NaCsPlot.maybe_save("$(prefix)_ratio")

# figure(figsize=[14.4, 4.5])
# plot(data[:, 1] ./ 1000, abs.(data[:, 2] ./ 2π / 1000), "C0", label="\$\\Omega_{R}\$")
# plot(data2[:, 1] ./ 1000, abs.(data2[:, 4] ./ 2π / 1000), "r", linewidth=4, alpha=0.7)
# plot(data[:, 1] ./ 1000, abs.(data[:, 4] ./ 2π / 1000), "C1", label="\$\\Gamma_{s}\$")
# ylabel("\$2\\pi\\cdot \\mathrm{kHz}\$")
# xlabel("Raman Single-Photon Frequency (THz)")
# yscale("log")
# ylim([0.12, 400])
# yticks([1, 10, 100], ["1", "10", "100"])
# grid()
# legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, 0.95), frameon=false)
# gca().yaxis.set_label_coords(-0.034, 0.55)
# xlim([287, 340])
# NaCsPlot.maybe_save("$(prefix)_full")

NaCsPlot.maybe_show()
