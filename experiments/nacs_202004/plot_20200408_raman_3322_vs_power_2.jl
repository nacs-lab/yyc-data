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

powers = [21, 18, 12, 0]
freqs = [772.1286, 771.924, 771.323, 770.202123]
freq_uncs = [0.0032, 0.033, 0.050, 0.000074]

model_pwr1(x, p) = p[1] .* (x)
model_pwr2(x, p) = p[1] .* (x).^2
model_pwr(x, p) = p[1] .* x.^p[2]
model_lin_off(x, p) = p[1] .+ p[2] .* x
model_sqr_off(x, p) = p[1] .+ p[2] .* x .+ p[3] .* x.^2

fit_freq = fit_data(model_lin_off, powers, freqs, freq_uncs, [770.2, 0.8])
fit_freq2 = fit_data(model_sqr_off, powers, freqs, freq_uncs, [770.2, 0.8, 0])

const prefix = joinpath(@__DIR__, "imgs", "data_20200408_raman_3322_vs_power_2")

figure()
errorbar(powers, freqs, freq_uncs, fmt="C0.")
plot(fit_freq.plotx, fit_freq.ploty, "C0", label="Linear")
plot(fit_freq2.plotx, fit_freq2.ploty, "C2", label="Quadratic")
text(0, 771.25, "\$f=f_0+a \\cdot P\$", color="C0", fontsize="small")
text(0, 771.55, ("\$f_0=$(fit_freq.uncs[1]) \\mathrm{MHz}\$\n" *
                 "\$a=$(fit_freq.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-1}\$"),
     color="C0", fontsize="small")
text(9.5, 770.83, "\$f=f_0+a \\cdot P+b \\cdot P^2\$", color="C2", fontsize="small")
text(7, 770.21, ("\$f_0=$(fit_freq2.uncs[1]) \\mathrm{MHz}\$\n" *
                 "\$a=$(fit_freq2.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-1}\$\n" *
                 "\$b=$(fit_freq2.uncs[3] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-2}\$"),
     color="C2", fontsize="small")
legend(fontsize="small", loc="upper left", ncol=2)
grid()
title("Raman Resonance (288668.35 GHz)")
xlabel("Tweezer Power (mW)")
ylabel("Raman Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)_f")

figure()
errorbar(powers, (freqs .- model_lin_off(powers, fit_freq.param)) .* 1000,
         freq_uncs .* 1000, fmt="C0.-", label="Linear")
errorbar(powers, (freqs .- model_sqr_off(powers, fit_freq2.param)) .* 1000,
         freq_uncs .* 1000, fmt="C2.-", label="Quadratic")
legend(fontsize="small")
grid()
title("Frequency Fit Residue")
xlabel("Tweezer Power (mW)")
ylabel("Residue (kHz)")
NaCsPlot.maybe_save("$(prefix)_fres")

NaCsPlot.maybe_show()
