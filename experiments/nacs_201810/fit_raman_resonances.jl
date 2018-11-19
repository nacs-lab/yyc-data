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
     "\$a=$(uncs[2])\\ \$MHz/mW\n" *
     "\$f_{PA0}=$(uncs[3])\\ \$GHz", fontsize="small")
text(370, 304.5, "\$f_{Raman}=f_{Raman0}-\\dfrac{a\\cdot P}{f_{PA} - f_{PA0}}\$")
grid()
xlim([350, 750])
ylim([293, 307])
title("Light shift (10mW)")
xlabel("288XXX (GHz)")
ylabel("Raman resonance (MHz)")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()

# exit()

# const times_a = Float64[0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]
# const ntimes_a = length(times_a)
# const freqs_a = 298.126:0.0005:298.134

# const spec_a = (1.0:(ntimes_a * length(freqs_a)), 0.0:0.0)
# const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)

# const time_datas =
#     [[split_nacs_a[2];
#       NaCsData.map_params((i, v)->times_a[i],
#                           split_nacs_a[1][((i - 1) * ntimes_a + 1):((i - 1) * ntimes_a + ntimes_a)])]
#      for i in 1:length(freqs_a)]
# const freq_datas =
#     [NaCsData.map_params((i, v)->freqs_a[i], split_nacs_a[1][i:ntimes_a:end])
#      for i in 1:ntimes_a]

# const prefix = joinpath(@__DIR__, "imgs", "data_20181114_raman_2d")

# figure()
# for data in time_datas
#     NaCsPlot.plot_survival_data(data, fmt=".-")
# end
# grid()
# ylim([0, 1])
# title("Raman time")
# xlabel("Time (\$ms\$)")
# ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)_time")

# figure()
# for data in freq_datas
#     NaCsPlot.plot_survival_data(data, fmt=".-")
# end
# grid()
# ylim([0, 1])
# title("Raman spectrum")
# xlabel("Detuning (\$MHz\$)")
# ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)_freq")

# NaCsPlot.maybe_show()
