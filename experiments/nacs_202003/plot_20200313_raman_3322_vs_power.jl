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

powers = [15, 12, 9, 6, 3]
freqs = [770.59511, 770.52272, 770.44752, 770.3699, 770.2875]
freq_uncs = [0.00030, 0.00079, 0.00091, 0.0010, 0.0014]

model_pwr1(x, p) = p[1] .* (x)
model_pwr2(x, p) = p[1] .* (x).^2
model_pwr(x, p) = p[1] .* x.^p[2]
model_lin_off(x, p) = p[1] .+ p[2] .* x
model_sqr_off(x, p) = p[1] .+ p[2] .* x .+ p[3] .* x.^2

fit_freq = fit_data(model_lin_off, powers, freqs, freq_uncs, [770.2, 25])
fit_freq2 = fit_data(model_sqr_off, powers, freqs, freq_uncs, [770.2, 25, 0])

const prefix = joinpath(@__DIR__, "imgs", "data_20200313_raman_3322_vs_power")

figure()
errorbar(powers, freqs, freq_uncs, fmt="C0.")
plot(fit_freq.plotx, fit_freq.ploty, "C0", label="Linear")
plot(fit_freq2.plotx, fit_freq2.ploty, "C2", label="Quadratic")
text(2.2, 770.46, "\$f=f_0+a \\cdot P\$", color="C0", fontsize="small")
text(2.2, 770.50, ("\$f_0=$(fit_freq.uncs[1]) \\mathrm{MHz}\$\n" *
                   "\$a=$(fit_freq.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-1}\$"),
     color="C0", fontsize="small")
text(9, 770.38, "\$f=f_0+a \\cdot P+b \\cdot P^2\$", color="C2", fontsize="small")
text(7, 770.27, ("\$f_0=$(fit_freq2.uncs[1]) \\mathrm{MHz}\$\n" *
                 "\$a=$(fit_freq2.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-1}\$\n" *
                 "\$b=$(fit_freq2.uncs[3] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-2}\$"),
     color="C2", fontsize="small")
legend(fontsize="small", loc="upper left", ncol=2)
grid()
title("Raman Resonance (288560 GHz)")
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
