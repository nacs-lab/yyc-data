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

# 1: 288503 GHz
# 2: 288560 GHz

const pwrs = [15, 6, 3]
const Γ_atom1 = [32.5, 2.77, 0.573]
const unc_Γ_atom1 = [3.1, 0.27, 0.078]
const Γ_atom2 = [24.5, 2.40, 0.572]
const unc_Γ_atom2 = [2.3, 0.20, 0.063]
const Γ_mol1 = [0.700, 0.227, 0.159]
const unc_Γ_mol1 = [0.045, 0.020, 0.023]
const Γ_mol2 = [1.30, 0.395, 0.170]
const unc_Γ_mol2 = [0.10, 0.041, 0.010]

const Γ_mol_diff = Γ_mol2 .- Γ_mol1
const unc_Γ_mol_diff = sqrt.(unc_Γ_mol1.^2 .+ unc_Γ_mol2.^2)

function gen_model_power(pwr)
    return (x, p)->x.^pwr .* p[1]
end

model_lin0(x, p) = x .* p[1]
model_lin1(x, p) = x .* p[1] .+ p[2]
model_quad0(x, p) = x.^2 .* p[1]

fit1_Γ_atom1 = fit_data(gen_model_power(1.58), pwrs, Γ_atom1, unc_Γ_atom1, [0.2], plot_lo=0)
fit1_Γ_atom2 = fit_data(gen_model_power(1.58), pwrs, Γ_atom2, unc_Γ_atom2, [0.2], plot_lo=0)
fit2_Γ_atom1 = fit_data(gen_model_power(2.58), pwrs, Γ_atom1, unc_Γ_atom1, [0.03], plot_lo=0)
fit2_Γ_atom2 = fit_data(gen_model_power(2.58), pwrs, Γ_atom2, unc_Γ_atom2, [0.03], plot_lo=0)

fit1_Γ_mol1 = fit_data(model_lin0, pwrs, Γ_mol1, unc_Γ_mol1, [0.05], plot_lo=0)
fit1_Γ_mol2 = fit_data(model_lin0, pwrs, Γ_mol2, unc_Γ_mol2, [0.05], plot_lo=0)
fit2_Γ_mol1 = fit_data(model_quad0, pwrs, Γ_mol1, unc_Γ_mol1, [0.01], plot_lo=0)
fit2_Γ_mol2 = fit_data(model_quad0, pwrs, Γ_mol2, unc_Γ_mol2, [0.01], plot_lo=0)

fit1_Γ_mol_diff = fit_data(model_lin0, pwrs, Γ_mol_diff, unc_Γ_mol_diff, [0.01], plot_lo=0)
fit2_Γ_mol_diff = fit_data(model_quad0, pwrs, Γ_mol_diff, unc_Γ_mol_diff, [0.01], plot_lo=0)

@show fit2_Γ_atom1.uncs # 2.93(17)\times10^{-2}
@show fit2_Γ_atom2.uncs # 2.45(13)\times10^{-2}
@show fit1_Γ_mol1.uncs # 4.35(21)\times10^{-2}
@show fit1_Γ_mol2.uncs # 6.32(27)\times10^{-2}

const prefix = joinpath(@__DIR__, "../figures/raman_transfer_scatter_scaling")

# figure()
# ax = gca()
# errorbar([], [], [], fmt="C0o-", label="288503 GHz")
# errorbar([], [], [], fmt="C1s-", label="288560 GHz")
# errorbar(pwrs, Γ_atom1, unc_Γ_atom1, fmt="C0o")
# plot(fit1_Γ_atom1.plotx, fit1_Γ_atom1.ploty, "C0--")
# plot(fit2_Γ_atom1.plotx, fit2_Γ_atom1.ploty, "C0")
# errorbar(pwrs, Γ_atom2, unc_Γ_atom2, fmt="C1s")
# plot(fit1_Γ_atom2.plotx, fit1_Γ_atom2.ploty, "C1--")
# plot(fit2_Γ_atom2.plotx, fit2_Γ_atom2.ploty, "C1")
# xlim([2.4, 17])
# ylim([0.4, 40])
# legend(loc="lower right", fontsize="x-small")
# xscale("log")
# yscale("log")
# grid()
# ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.set_xticks(3:17, minor=true)
# xticks([3, 6, 12], ["3", "6", "12"])
# yticks([1, 3, 10, 30], ["1", "", "10", ""])
# xlabel("Tweezer Power (mW)")
# ylabel("\$\\Gamma_{\\mathrm{a}}~(\\mathrm{2\\pi\\times Hz})\$")
# NaCsPlot.maybe_save("$(prefix)_atom")

# figure()
# ax = gca()
# errorbar([], [], [], fmt="C1s-", label="288560 GHz")
# errorbar([], [], [], fmt="C0o-", label="288503 GHz")
# errorbar(pwrs, Γ_mol1, unc_Γ_mol1, fmt="C0o")
# # plot(fit2_Γ_mol1.plotx, fit2_Γ_mol1.ploty, "C0--")
# plot(fit1_Γ_mol1.plotx, fit1_Γ_mol1.ploty, "C0")
# errorbar(pwrs, Γ_mol2, unc_Γ_mol2, fmt="C1s")
# # plot(fit2_Γ_mol2.plotx, fit2_Γ_mol2.ploty, "C1--")
# plot(fit1_Γ_mol2.plotx, fit1_Γ_mol2.ploty, "C1")
# xlim([2.4, 17])
# ylim([0.1, 1.6])
# legend(loc="lower right", fontsize="x-small")
# xscale("log")
# yscale("log")
# grid()
# ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.set_xticks(3:17, minor=true)
# xticks([3, 6, 12], ["3", "6", "12"])
# yticks([0.2, 0.4, 0.8, 1.6], ["0.2", "0.4", "0.8", "1.6"])
# xlabel("Tweezer Power (mW)")
# ylabel("\$\\Gamma_{\\mathrm{m}}~(\\mathrm{2\\pi\\times kHz})\$")
# NaCsPlot.maybe_save("$(prefix)_mol")

const diffscale = (1 / (711 - 560)^2) / (1 / (711 - 560)^2 - 1 / (711 - 503)^2)

figure()
ax = gca()
errorbar(pwrs, Γ_mol_diff .* diffscale, unc_Γ_mol_diff .* diffscale, fmt="C0o")
plot(fit1_Γ_mol_diff.plotx, fit1_Γ_mol_diff.ploty .* diffscale, "C0--", label="Linear")
plot(fit2_Γ_mol_diff.plotx, fit2_Γ_mol_diff.ploty .* diffscale, "C0", label="Quadratic")
xlim([2.4, 17])
ylim([0.016, 1.6])
xscale("log")
yscale("log")
legend(loc="upper left", fontsize="small")
grid()
text(3.1, 0.031,
     ("\$\\Gamma^{\\mathrm{\\Delta dep}}_{\\mathrm{scatter}}=P^2\\cdot2\\pi\\times$(fit2_Γ_mol_diff.uncs[1] .* 1000 * diffscale)" *
      "~\\mathrm{Hz/mW^2}\$"), fontsize="small", color="C0")
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks(3:17, minor=true)
xticks([3, 6, 12], ["3", "6", "12"])
yticks([0.03, 0.1, 0.3, 1], ["", "0.1", "", "1"])
xlabel("Tweezer Power (mW)")
ylabel("\$\\Gamma^{\\mathrm{\\Delta dep}}_{\\mathrm{scatter}}~(\\mathrm{2\\pi\\times kHz})\$")
gca().yaxis.set_label_coords(-0.09, 0.55)
NaCsPlot.maybe_save("$(prefix)_mol_diff")

NaCsPlot.maybe_show()
