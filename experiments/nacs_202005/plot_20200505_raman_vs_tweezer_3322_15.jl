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

const tweezer = [668.35, 605, 565, 560, 555, 510, 503, 492]
const freq = [1623, 760.9, 616.56, 595.80, 577.44, 501.17, 483.50, 475.71]
const freq_unc = [30, 1.0, 0.21, 0.47, 0.13, 0.83, 0.12, 0.10]

const prefix = joinpath(@__DIR__, "imgs", "data_20200505_ramann_vs_tweezer_3322_15")

function tweezer_model(x, p)
    return p[1] .- p[2] ./ (x .- p[3])
end

function gen_tweezer_model_res(f)
    return function(x, p)
        return p[1] .- p[2] ./ (x .- f)
    end
end

fit_freq = fit_data(tweezer_model, tweezer, freq, freq_unc, [250.0, 60000, 705],
                    plot_lo=300, plot_hi=670)
fit_freq2 = fit_data(gen_tweezer_model_res(705),
                     tweezer[2:end], freq[2:end], freq_unc[2:end], [250.0, 60000],
                     plot_lo=300, plot_hi=670)
@show fit_freq.uncs
@show fit_freq2.uncs

figure()
errorbar(tweezer, freq, freq_unc, fmt="C0s")
plot(fit_freq.plotx, fit_freq.ploty, "C0", label="With 668.35 GHz point")
plot(fit_freq2.plotx, fit_freq2.ploty, "C1", label="\$f_{res}=705\\ GHz\$")
axhline(202, color="C3", ls="--")
legend(fontsize="small")
grid()
xlabel("Tweezer Frequency (288XXX GHz)")
ylabel("Raman Resonance (770XXX kHz)")
title("Raman Resonance (15 mW)")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
