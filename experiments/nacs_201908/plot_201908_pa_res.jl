#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
using MAT

# 15mW, 111(11) MHz, 203(41) MHz
# 15mW, 1107(61) MHz, 220(29) MHz
# 10mW, 98.7(62) MHz, 167(32) MHz
# 5mW, 63.9(32) MHz, 50(15) MHz
# 5mW, 64.4(50) MHz, 64(22) MHz
# 1mW, 38.6(19) MHz, 27.3(97) MHz

# Γ = 13.9(11) MHz/mW * tweezer_power

const pwrs = [15, 10, 5, 1]
const freqs = [110.8, 98.7, 64.0, 38.6]
const freqs_s = [5.3, 6.2, 2.6, 1.9]
const Γs = [214, 167, 54, 27.3]
const Γs_s = [23, 32, 12, 9.7]

model_lin0(x, p) = x .* p[1]
model_lin1(x, p) = x .* p[1] .+ p[2]

fit_freq = fit_data(model_lin1, pwrs, freqs, freqs_s, [6.0, 30], plot_lo=0)
fit_Γ = fit_data(model_lin0, pwrs, Γs, Γs_s, [10.0], plot_lo=0)

const prefix = joinpath(@__DIR__, "imgs", "data_201908_pa_res")

figure()
errorbar(pwrs, freqs, freqs_s, fmt="C0.")
plot(fit_freq.plotx, fit_freq.ploty, "C0")
text(1, 100, "\$f = f_0 + b \\cdot P\$", color="C0", fontsize="large")
text(5, 32, ("\$f_0 = $(fit_freq.uncs[2]) \\mathrm{MHz}\$\n" *
             "\$b = $(fit_freq.uncs[1]) \\mathrm{MHz\\cdot mW^{-1}}\$"), color="C0")
xlim([0, 16])
grid()
xlabel("Tweezer Power (mW)")
ylabel("306496XXX MHz")
title("PA resonance")
NaCsPlot.maybe_save("$(prefix)_freq")

figure()
errorbar(pwrs, Γs, Γs_s, fmt="C0.")
plot(fit_Γ.plotx, fit_Γ.ploty, "C0")
text(1, 200, "\$\\Gamma = b \\cdot P\$", color="C0", fontsize="large")
text(5, 10, "\$b = $(fit_Γ.uncs[1]) \\mathrm{MHz\\cdot mW^{-1}}\$", color="C0")
xlim([0, 16])
ylim([0, 260])
grid()
xlabel("Tweezer Power (mW)")
ylabel("Linewidth (MHz)")
title("PA linewidth")
NaCsPlot.maybe_save("$(prefix)_gamma")

NaCsPlot.maybe_show()
