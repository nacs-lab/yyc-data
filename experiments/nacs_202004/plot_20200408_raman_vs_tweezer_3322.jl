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

const tweezer = [668.35, 565, 560, 555]
const freq = [488, 291.891, 287.282, 283.665]
const freq_unc = [10, 0.038, 0.015, 0.017]

const prefix = joinpath(@__DIR__, "imgs", "data_20200408_raman_vs_tweezer_3322")

function tweezer_model(x, p)
    return p[1] .- p[2] ./ (x .- p[3])
end

function gen_tweezer_model_res(f)
    return function(x, p)
        return p[1] .- p[2] ./ (x .- f)
    end
end

function gen_tweezer_model_res2(f1, f2)
    return function(x, p)
        return p[1] .- p[2] ./ (x .- f1) .- p[3] ./ (x .- f2)
    end
end

fit_freq = fit_data(tweezer_model, tweezer, freq, freq_unc, [250.0, 60000, 705],
                    plot_lo=300, plot_hi=670)
fit_freq2 = fit_data(gen_tweezer_model_res(705),
                     tweezer[2:end], freq[2:end], freq_unc[2:end], [250.0, 60000],
                     plot_lo=300, plot_hi=670)
fit_freq3 = fit_data(gen_tweezer_model_res2(705, -15),
                     tweezer, freq, freq_unc, [250.0, 60000, 10000],
                     plot_lo=300, plot_hi=670)
@show fit_freq.uncs
@show fit_freq2.uncs
@show fit_freq3.uncs

figure()
errorbar(tweezer, freq, freq_unc, fmt="C0s")
plot(fit_freq.plotx, fit_freq.ploty, "C0", label="With 668.35 GHz point")
plot(fit_freq2.plotx, fit_freq2.ploty, "C1", label="\$f_{res}=705\\ GHz\$")
plot(fit_freq3.plotx, fit_freq3.ploty, "C2", label="\$f_{res}=705, -15\\ GHz\$")
legend(fontsize="small")
grid()
xlabel("Tweezer Frequency (288XXX GHz)")
ylabel("Raman Resonance (770XXX kHz)")
title("Raman Resonance (3 mW)")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
