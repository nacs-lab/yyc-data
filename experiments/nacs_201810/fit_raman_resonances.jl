#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
using DelimitedFiles

const iname = joinpath(@__DIR__, "data", "raman_resonance.csv")
const rdata = readdlm(iname, ',', Float64, skipstart=1)
const single_freqs = rdata[:, 1] .+ rdata[:, 5] / 1000
const powers = rdata[:, 2]
const res_freqs = rdata[:, 3]
const res_freqs_unc = rdata[:, 4] .* 2

function resonance(sfs, pwrs, p)
    return p[1] .+ p[2] .* pwrs ./ (p[3] .- sfs)
end

function model(x, p)
    sfs = single_freqs[x]
    pwrs = powers[x]
    return resonance(sfs, pwrs, p)
end

fit = curve_fit(model, 1:length(single_freqs), res_freqs, res_freqs_unc.^-(2/3),
                [298.0, 7.5, 694.5])
uncs = Unc.(fit.param, estimate_errors(fit))
@show uncs

const prefix = joinpath(@__DIR__, "imgs", "fit_raman_resonances")

plotf1 = linspace(350, fit.param[3] - 0.1, 1000)
plotf2 = linspace(fit.param[3] + 0.1, 750, 1000)

figure()
plot(plotf1, resonance(plotf1, 10, fit.param), "C1")
plot(plotf2, resonance(plotf2, 10, fit.param), "C1")
errorbar(single_freqs[powers .== 10], res_freqs[powers .== 10], res_freqs_unc[powers .== 10],
         fmt="C0.")
text(380, 299.5, "\$f_{Raman0}=$(uncs[1])\\ \$MHz\n" *
     "\$a=$(uncs[2])\\ \$MHz\$\\cdot\$GHz/mW\n" *
     "\$f_{PA0}=$(uncs[3])\\ \$GHz", fontsize="small")
text(370, 304.5, "\$f_{Raman}=f_{Raman0}-\\dfrac{a\\cdot P}{f_{PA} - f_{PA0}}\$")
grid()
xlim([350, 750])
ylim([293, 307])
title("Light shift (10mW)")
xlabel("288XXX (GHz)")
ylabel("Raman resonance (MHz)")
NaCsPlot.maybe_save("$(prefix)")

figure()
function plot_power(p, plotfs, color)
    first = true
    for plotf in plotfs
        if first
            plot(plotf, resonance(plotf, p, fit.param), color, label="$p mW")
            first = false
        else
            plot(plotf, resonance(plotf, p, fit.param), color)
        end
    end
    errorbar(single_freqs[powers .== p], res_freqs[powers .== p], res_freqs_unc[powers .== p],
             fmt="$color.")
end

plot_power(2, (linspace(fit.param[3] + 0.1, 750, 1000),), "C0")
plot_power(4, (linspace(fit.param[3] + 0.1, 750, 1000),), "C1")
plot_power(10, (plotf1, plotf2), "C2")
plot_power(15, (linspace(540, 625, 1000),), "C3")
plot_power(20, (linspace(540, 625, 1000),), "C4")
text(380, 300, "\$f_{Raman0}=$(uncs[1])\\ \$MHz\n" *
     "\$a=$(uncs[2])\\ \$MHz\$\\cdot\$GHz/mW\n" *
     "\$f_{PA0}=$(uncs[3])\\ \$GHz", fontsize="small")
text(370, 304.5, "\$f_{Raman}=f_{Raman0}-\\dfrac{a\\cdot P}{f_{PA} - f_{PA0}}\$")
legend(ncol=2, fontsize="small", loc="lower left", labelspacing=0.4)
grid()
xlim([350, 750])
ylim([293, 307])
title("Light shift")
xlabel("288XXX (GHz)")
ylabel("Raman resonance (MHz)")
NaCsPlot.maybe_save("$(prefix)_all")

NaCsPlot.maybe_show()
