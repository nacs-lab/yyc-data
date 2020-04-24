#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200421_223616.mat",
                "data_20200422_093717.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [(484.00 .+ [-30; -4.0:0.8:4.0; 30], # 15 mW, 0.20 ms
                433.00 .+ [-20; -2.5:0.5:2.5; 20], # 12 mW, 0.25 ms
                379.75 .+ [-10; -1.5:0.3:1.5; 10], #  9 mW, 0.40 ms
                324.5 .+ [-8; -1.25:0.25:1.25; 8], #  6 mW, 0.80 ms
                265.00 .+ [-5; -1.0:0.2:1.0; 5], # 3 mW, 2.0 ms
                ),
               (484.00 .+ [-30; -4.0:0.8:4.0; 30], # 15 mW, 0.20 ms
                433.00 .+ [-20; -2.5:0.5:2.5; 20], # 12 mW, 0.25 ms
                379.75 .+ [-10; -1.5:0.3:1.5; 10], #  9 mW, 0.40 ms
                324.5 .+ [-8; -1.25:0.25:1.25; 8], #  6 mW, 0.80 ms
                265.00 .+ [-5; -1.0:0.2:1.0; 5], # 3 mW, 2.0 ms
                )
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const data_nacs_15 = [datas_nacs[1][1]; datas_nacs[2][1]]
const data_nacs_12 = [datas_nacs[1][2]; datas_nacs[2][2]]
const data_nacs_9 = [datas_nacs[1][3]; datas_nacs[2][3]]
const data_nacs_6 = [datas_nacs[1][4]; datas_nacs[2][4]]
const data_nacs_3 = [datas_nacs[1][5]; datas_nacs[2][5]]

const prefix = joinpath(@__DIR__, "imgs", "data_20200421_223616_raman_det_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
model_pwr1(x, p) = p[1] .* (x)
model_pwr2(x, p) = p[1] .* (x).^2
model_pwr(x, p) = p[1] .* x.^p[2]
model_lin_off(x, p) = p[1] .+ p[2] .* x
model_sqr_off(x, p) = p[1] .+ p[2] .* x .+ p[3] .* x.^2

fit15 = fit_survival(model_lorentzian, data_nacs_15, [0.4, 0.3, 484.0, 8])
fit12 = fit_survival(model_lorentzian, data_nacs_12, [0.4, 0.3, 433.0, 5])
fit9 = fit_survival(model_lorentzian, data_nacs_9, [0.4, 0.3, 379.75, 3])
fit6 = fit_survival(model_lorentzian, data_nacs_6, [0.3, 0.25, 324.5, 2.5])
fit3 = fit_survival(model_lorentzian, data_nacs_3, [0.4, 0.3, 265.0, 2])

powers = [15, 12, 9, 6, 3]
freqs = [770 + fit.param[3] / 1000 for fit in [fit15, fit12, fit9, fit6, fit3]]
freq_uncs = [fit.unc[3] / 1000 for fit in [fit15, fit12, fit9, fit6, fit3]]

fit_freq = fit_data(model_lin_off, powers, freqs, freq_uncs, [770.2, 25])
fit_freq2 = fit_data(model_sqr_off, powers, freqs, freq_uncs, [770.2, 25, 0])

@show fit15.uncs
@show fit12.uncs
@show fit9.uncs
@show fit6.uncs
@show fit3.uncs

figure(figsize=[12.6, 16.8])

subplot(3, 2, 1)
NaCsPlot.plot_survival_data(data_nacs_15, fmt="C0.")
plot(fit15.plotx, fit15.ploty, "C0")
text(460, 0.332, "\$f=$(770 + fit15.uncs[3] / 1000)\\ MHz\$")
title("288503 GHz, 15 mW, 0.19 ms")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")

subplot(3, 2, 2)
NaCsPlot.plot_survival_data(data_nacs_12, fmt="C0.")
plot(fit12.plotx, fit12.ploty, "C0")
text(415, 0.280, "\$f=$(770 + fit12.uncs[3] / 1000)\\ MHz\$")
title("288503 GHz, 12 mW, 0.24 ms")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")

subplot(3, 2, 3)
NaCsPlot.plot_survival_data(data_nacs_9, fmt="C0.")
plot(fit9.plotx, fit9.ploty, "C0")
text(372.5, 0.270, "\$f=$(770 + fit9.uncs[3] / 1000)\\ MHz\$")
title("288503 GHz, 9 mW, 0.39 ms")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")

subplot(3, 2, 4)
NaCsPlot.plot_survival_data(data_nacs_6, fmt="C0.")
plot(fit6.plotx, fit6.ploty, "C0")
text(319, 0.299, "\$f=$(770 + fit6.uncs[3] / 1000)\\ MHz\$")
title("288503 GHz, 6 mW, 0.79 ms")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")

subplot(3, 2, 5)
NaCsPlot.plot_survival_data(data_nacs_3, fmt="C0.")
plot(fit3.plotx, fit3.ploty, "C0")
text(260.5, 0.328, "\$f=$(770 + fit3.uncs[3] / 1000)\\ MHz\$")
title("288503 GHz, 3 mW, 1.99 ms")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

figure()
errorbar(powers, freqs, freq_uncs, fmt="C0.")
plot(fit_freq.plotx, fit_freq.ploty, "C0", label="Linear")
plot(fit_freq2.plotx, fit_freq2.ploty, "C2", label="Quadratic")
text(2.05, 770.355, "\$f=f_0+a \\cdot P\$", color="C0", fontsize="small")
text(2.05, 770.415, ("\$f_0=$(fit_freq.uncs[1]) \\mathrm{MHz}\$\n" *
                     "\$a=$(fit_freq.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-1}\$"),
     color="C0", fontsize="small")
text(9, 770.352, "\$f=f_0+a \\cdot P+b \\cdot P^2\$", color="C2", fontsize="small")
text(6.5, 770.25, ("\$f_0=$(fit_freq2.uncs[1]) \\mathrm{MHz}\$\n" *
                   "\$a=$(fit_freq2.uncs[2] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-1}\$\n" *
                   "\$b=$(fit_freq2.uncs[3] * 1000) \\mathrm{kHz}\\cdot\\mathrm{mW}^{-2}\$"),
     color="C2", fontsize="small")
legend(fontsize="small", loc="upper left", ncol=2)
grid()
title("Raman Resonance (288503 GHz)")
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
