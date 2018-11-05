#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181103_165031.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20181104_020816.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const iname_c = joinpath(@__DIR__, "data", "data_20181104_134413.mat")
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)
const iname_d = joinpath(@__DIR__, "data", "data_20181105_002615.mat")
const params_d, logicals_d = NaCsData.load_striped_mat(iname_d)
const iname_e = joinpath(@__DIR__, "data", "data_20181105_075855.mat")
const params_e, logicals_e = NaCsData.load_striped_mat(iname_e)

data_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2), (3, 4)))
data_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((1, 2), (3, 4)))
data_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((1, 2), (3, 4)))
data_d = NaCsData.select_count(params_d, logicals_d, NaCsData.select_single((1, 2), (3, 4)))
data_e = NaCsData.select_count(params_e, logicals_e, NaCsData.select_single((1, 2), (3, 4)))

const spec_a = (81:0.06:90) .* 2 .- 60
const spec_b = (79:0.06:88) .* 2 .- 60
const spec_c = ((81.9:0.06:85.5) .* 2 .- 60,
                (82.1:0.06:86) .* 2 .- 60,
                (83:0.06:89) .* 2 .- 60)
const spec_d = ((81.1:0.05:84.2) .* 2 .- 60,
                (82.62:0.06:85.02) .* 2 .- 60,
                (83:0.06:85.76) .* 2 .- 60,
                (83.66:0.06:87.5) .* 2 .- 60)
const spec_e = (81.6:0.05:84.2) .* 2 .- 60

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)
const split_c = NaCsData.split_data(data_c, spec_c)
const split_d = NaCsData.split_data(data_d, spec_d)
const split_e = NaCsData.split_data(data_e, spec_e)

data_100 = split_a
data_60 = split_b
data_73 = [split_c[1]; split_d[2]]
data_86 = [split_c[2]; split_d[3]]
data_120 = [split_c[3]; split_d[4]]
data_40 = [split_d[1]; split_e]

const freqs = [4 103.63 0.05
               4 105.70 0.05
               4 107.70 0.1
               6 105.32 0.06
               6 107.18 0.06
               6 108.164 0.06
               7.3 106.2 0.06
               7.3 108 0.06
               7.3 109 0.06
               8.6 106.74 0.06
               8.6 108.4 0.06
               8.6 110.26 0.10
               10 107.34 0.06
               10 109.02 0.06
               10 111.72 0.10
               12 107.8 0.10
               12 109.864 0.06
               12 113.92 0.10]
const powers = [4, 6, 7.3, 8.6, 10, 12]
function elevels(power, ps)
    A = [ps[1] + power * ps[4] ps[7] ps[8]
         ps[7] ps[2] + power * ps[5] ps[9]
         ps[8] ps[9] ps[3] + power * ps[6]]
    return eigvals(A)
end

function model(x, ps)
    Es = Float64[]
    for power in powers
        append!(Es, elevels(power, ps))
    end
    return Es[x]
end
fit = curve_fit(model, 1:size(freqs, 1), freqs[:, 2], [99.7, 103.1, 106.9,
                                                       1.16, 0.55, 0.11,
                                                       -0.5, 0.9, 0.6])
uncs = Unc.(fit.param, estimate_errors(fit))
@show uncs

function eval_model(powers, ps)
    E1s = Float64[]
    E2s = Float64[]
    E3s = Float64[]
    for power in powers
        E = elevels(power, ps)
        append!(E1s, E[1])
        append!(E2s, E[2])
        append!(E3s, E[3])
    end
    return powers, E1s, E2s, E3s
end

const prefix = joinpath(@__DIR__, "imgs", "data_20181103_raman_n2")

figure()
NaCsPlot.plot_survival_data(data_60, fmt="C0.-", label="6mW")
NaCsPlot.plot_survival_data(data_100, fmt="C3.-", label="10mW")
grid()
legend()
ylim([0, 1])
title("N=2 Raman")
xlabel("Raman Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_6_10")

figure()
NaCsPlot.plot_survival_data(data_40, fmt="C0.-", label="4mW")
NaCsPlot.plot_survival_data(data_60, fmt="C1.-", label="6mW")
NaCsPlot.plot_survival_data(data_73, fmt="C2.-", label="7.3mW")
grid()
legend()
ylim([0, 1])
title("N=2 Raman")
xlabel("Raman Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_40-73")

figure()
NaCsPlot.plot_survival_data(data_86, fmt="C0.-", label="8.6mW")
NaCsPlot.plot_survival_data(data_100, fmt="C1.-", label="10mW")
NaCsPlot.plot_survival_data(data_120, fmt="C2.-", label="12mW")
grid()
legend()
ylim([0, 1])
title("N=2 Raman")
xlabel("Raman Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_86-120")

plotx, plote1, plote2, plote3 = eval_model(linspace(3.5, 13, 1000), fit.param)
uparams = copy(fit.param)
uparams[7:9] = 0
plotx, plotue1, plotue2, plotue3 = eval_model(linspace(3.5, 13, 1000), uparams)

@show elevels(0, fit.param)

figure()
plot(plotx, plote1, "C1", label="Fit")
plot(plotx, plote2, "C1")
plot(plotx, plote3, "C1")
plot(plotx, plotue1, "C2", ls="dotted", label="Uncoupled")
plot(plotx, plotue2, "C2", ls="dotted")
plot(plotx, plotue3, "C2", ls="dotted")
errorbar(freqs[:, 1], freqs[:, 2], freqs[:, 3], fmt="C0o", label="Measured")
grid()
legend(fontsize=18)
title("N=2 Raman Freq vs Power")
xlabel("Power (mW)")
ylabel("Raman Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)_powers")

params2 = copy(fit.param)
params2[2] -= 2.8 / 6 * 8.85 * 0.4
params2[3] -= 2.8 / 6 * 8.85 * 0.4 * 2
plot2x, plot2e1, plot2e2, plot2e3 = eval_model(linspace(0, 10, 1000), params2)

figure()
plot(plot2x, plot2e1, "C0")
plot(plot2x, plot2e2, "C0")
plot(plot2x, plot2e3, "C0")
# plot(plotx, plote1, "C1")
# plot(plotx, plote2, "C1")
# plot(plotx, plote3, "C1")
grid()
xlim([0, 11])
title("Prediction for 60% B field")
xlabel("Power (mW)")
ylabel("Raman Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)_b60")

NaCsPlot.maybe_show()
