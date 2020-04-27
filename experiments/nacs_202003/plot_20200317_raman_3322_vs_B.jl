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

Bs = [100, 95, 80, 60, 40]
freqs = [770.59511, 770.5766, 770.521, 770.4479, 770.3713]
freq_uncs = [0.00030, 0.0010, 0.004, 0.0014, 0.0017]

model_pwr1(x, p) = p[1] .* (x)
model_pwr2(x, p) = p[1] .* (x).^2
model_pwr(x, p) = p[1] .* x.^p[2]
model_lin_off(x, p) = p[1] .+ p[2] .* x
model_sqr_off(x, p) = p[1] .+ p[2] .* x .+ p[3] .* x.^2

fit_freq = fit_data(model_lin_off, Bs, freqs, freq_uncs, [770.2, 0.003])
# fit_freq2 = fit_data(model_sqr_off, Bs, freqs, freq_uncs, [770.2, 25, 0])

const prefix = joinpath(@__DIR__, "imgs", "data_20200317_raman_3322_vs_B")

figure()
errorbar(Bs, freqs, freq_uncs, fmt="C0.")
plot(fit_freq.plotx, fit_freq.ploty, "C0", label="Linear")
# plot(fit_freq2.plotx, fit_freq2.ploty, "C2", label="Quadratic")
text(62, 770.41, "\$f=f_0+a \\cdot B\$", color="C0", fontsize="small")
text(55, 770.355, ("\$f_0=$(fit_freq.uncs[1]) \\mathrm{MHz}\$\n" *
                   "\$a=$(fit_freq.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{\\%}^{-1}\$"),
     color="C0", fontsize="small")
# text(9, 770.38, "\$f=f_0+a \\cdot B+b \\cdot B^2\$", color="C2", fontsize="small")
# text(7, 770.27, ("\$f_0=$(fit_freq2.uncs[1]) \\mathrm{MHz}\$\n" *
#                  "\$a=$(fit_freq2.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{\\%}^{-1}\$\n" *
#                  "\$b=$(fit_freq2.uncs[3] * 1000) \\mathrm{kHz}\\cdot\\mathrm{\\%}^{-2}\$"),
#      color="C2", fontsize="small")
# legend(fontsize="small", loc="upper left", ncol=2)
grid()
title("Raman Resonance 288560 GHz 15 mW")
xlabel("B field (%)")
ylabel("Raman Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)")

figure()
errorbar(Bs, (freqs .- model_lin_off(Bs, fit_freq.param)) .* 1000,
         freq_uncs .* 1000, fmt="C0.-", label="Linear")
# errorbar(Bs, (freqs .- model_sqr_off(Bs, fit_freq2.param)) .* 1000,
#          freq_uncs .* 1000, fmt="C2.-", label="Quadratic")
# legend(fontsize="small")
grid()
title("Frequency Fit Residue")
xlabel("B field (%)")
ylabel("Residue (kHz)")
NaCsPlot.maybe_save("$(prefix)_res")

NaCsPlot.maybe_show()
