#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using PyCall
using DataStructures
using LsqFit
using MAT

const Γ_mol_0 = [1.87, 0.528, 0.298]
const unc_Γ_mol_0 = [0.10, 0.045, 0.028]
const Γ_mol_1 = [1.30, 0.395, 0.170]
const unc_Γ_mol_1 = [0.10, 0.041, 0.010]
const Γ_mol_2 = [0.781]
const unc_Γ_mol_2 = [0.050]

const Ω_mol_0 = [4.023, 1.182, 0.488]
const unc_Ω_mol_0 = [0.059, 0.033, 0.015]
const Ω_mol_1 = [3.971, 1.173, 0.4912]
const unc_Ω_mol_1 = [0.089, 0.025, 0.0082]
const Ω_mol_2 = [3.282]
const unc_Ω_mol_2 = [0.042]

const prefix = joinpath(@__DIR__, "../figures/raman_transfer_gamma_vs_omega")

figure()
ax = gca()
errorbar(Ω_mol_0, Γ_mol_0, xerr=unc_Ω_mol_0, yerr=unc_Γ_mol_0, fmt="C3o-", label="0 filter")
errorbar(Ω_mol_1, Γ_mol_1, xerr=unc_Ω_mol_1, yerr=unc_Γ_mol_1, fmt="C0o-", label="1 filter")
# errorbar(Ω_mol_2, Γ_mol_2, xerr=unc_Ω_mol_2, yerr=unc_Γ_mol_2, fmt="C2o-", label="2 filter")
legend(loc="upper center", fontsize="x-small")
xscale("log")
yscale("log")
grid()
xlim([0.42, 4.4])
ylim([0.14, 2.2])
minorticks_off()
xticks([0.5, 1, 2, 4], ["0.5", "1", "2", "4"])
yticks([0.2, 0.4, 0.8, 1.6], ["0.2", "0.4", "0.8", "1.6"])
# ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
xlabel("\$\\Omega_{\\mathrm{Raman}}~(\\mathrm{2\\pi\\times kHz})\$")
ylabel("\$\\Gamma_{\\mathrm{scatter}}~(\\mathrm{2\\pi\\times kHz})\$")
NaCsPlot.maybe_save("$(prefix)_0-1")

figure()
ax = gca()
errorbar(Ω_mol_0, Γ_mol_0, xerr=unc_Ω_mol_0, yerr=unc_Γ_mol_0, fmt="C3o-", label="0 filter")
errorbar(Ω_mol_1, Γ_mol_1, xerr=unc_Ω_mol_1, yerr=unc_Γ_mol_1, fmt="C0o-", label="1 filter")
errorbar(Ω_mol_2, Γ_mol_2, xerr=unc_Ω_mol_2, yerr=unc_Γ_mol_2, fmt="C2o-", label="2 filter")
annotate("", xy=(Ω_mol_2[1] - 0.3 / 4, Γ_mol_2[1] - 0.3 / 4),
         xytext=(Ω_mol_2[1] - 0.3 * 1.25, Γ_mol_2[1] - 0.3 * 1.25),
         arrowprops=Dict(:width=>5, :headlength=>10, :headwidth=>10, :color=>"C2"))
legend(loc="upper center", fontsize="x-small")
xscale("log")
yscale("log")
grid()
xlim([0.42, 4.4])
ylim([0.14, 2.2])
minorticks_off()
xticks([0.5, 1, 2, 4], ["0.5", "1", "2", "4"])
yticks([0.2, 0.4, 0.8, 1.6], ["0.2", "0.4", "0.8", "1.6"])
# ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
xlabel("\$\\Omega_{\\mathrm{Raman}}~(\\mathrm{2\\pi\\times kHz})\$")
ylabel("\$\\Gamma_{\\mathrm{scatter}}~(\\mathrm{2\\pi\\times kHz})\$")
NaCsPlot.maybe_save("$(prefix)_0-2")

NaCsPlot.maybe_show()
