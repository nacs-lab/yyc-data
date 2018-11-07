#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181105_171851.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20181105_185338.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const iname_c = joinpath(@__DIR__, "data", "data_20181105_222500.mat")
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)
const iname_d = joinpath(@__DIR__, "data", "data_20181105_232836.mat")
const params_d, logicals_d = NaCsData.load_striped_mat(iname_d)
const iname_e = joinpath(@__DIR__, "data", "data_20181106_084221.mat")
const params_e, logicals_e = NaCsData.load_striped_mat(iname_e)

data_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2), (3, 4)))
data_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((1, 2), (3, 4)))
data_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((1, 2), (3, 4)))
data_d = NaCsData.select_count(params_d, logicals_d, NaCsData.select_single((1, 2), (3, 4)))
data_e = NaCsData.select_count(params_e, logicals_e, NaCsData.select_single((1, 2), (3, 4)))

const spec_a = (81.0:0.05:83.0) .* 2 .- 60 # 4
const spec_b = (80.5:0.04:82.5) .* 2 .- 60 # 2
const spec_c = (80.6:0.06:86.0) .* 2 .- 60 # 8
const spec_d = ((81.0:0.05:82.7) .* 2 .- 60, # 3
                (81.3:0.06:83.1) .* 2 .- 60, # 5
                (81.5:0.06:83.6) .* 2 .- 60, # 6
                (81.6:0.06:84.3) .* 2 .- 60, # 7
                (80.6:0.06:86.0) .* 2 .- 60) # 8
const spec_e = ((80.5:0.04:82.5) .* 2 .- 60, # 2
                (81.0:0.05:82.7) .* 2 .- 60, # 3
                (81.0:0.05:83.0) .* 2 .- 60, # 4
                (81.3:0.06:83.1) .* 2 .- 60, # 5
                (81.5:0.06:83.6) .* 2 .- 60, # 6
                (81.6:0.06:84.3) .* 2 .- 60, # 7
                (80.6:0.06:86.0) .* 2 .- 60) # 8

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)
const split_c = NaCsData.split_data(data_c, spec_c)
const split_d = NaCsData.split_data(data_d, spec_d)
const split_e = NaCsData.split_data(data_e, spec_e)

data_20 = split_b
data_30 = [split_d[1]; split_e[2]]
data_40 = split_a
data_50 = [split_d[2]; split_e[4]]
data_60 = [split_d[3]; split_e[5]]
data_70 = [split_d[4]; split_e[6]]
data_80 = [split_c; split_d[5]; split_e[7]]

const freqs = [2 101.480 0.05
               2 102.816 0.05
               2 104.500 0.12
               3 102.360 0.05
               3 103.600 0.05
               3 104.400 0.05
               4 103.050 0.05
               4 104.200 0.05
               4 104.860 0.05
               5 103.560 0.06
               5 104.688 0.06
               5 105.840 0.06
               6 103.924 0.06
               6 105.088 0.06
               6 106.960 0.06
               7 104.160 0.12
               7 105.600 0.06
               7 107.800 0.12
               8 104.320 0.12
               8 106.060 0.12
               8 108.800 0.15]
const powers = [2, 3, 4, 5, 6, 7, 8]
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
fit = curve_fit(model, 1:size(freqs, 1), freqs[:, 2], [99.8, 101.3, 103.9,
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

const prefix = joinpath(@__DIR__, "imgs", "data_20181105_raman_n2_b60")

figure()
NaCsPlot.plot_survival_data(data_20, fmt="C0.-", label="2mW")
NaCsPlot.plot_survival_data(data_30, fmt="C1.-", label="3mW")
NaCsPlot.plot_survival_data(data_40, fmt="C2.-", label="4mW")
grid()
legend()
ylim([0, 0.5])
title("N=2 Raman")
xlabel("Raman Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_2-4")

figure()
NaCsPlot.plot_survival_data(data_50, fmt="C0.-", label="5mW")
NaCsPlot.plot_survival_data(data_60, fmt="C1.-", label="6mW")
grid()
legend()
ylim([0, 0.5])
title("N=2 Raman")
xlabel("Raman Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_5-6")

figure()
NaCsPlot.plot_survival_data(data_70, fmt="C0.-", label="7mW")
NaCsPlot.plot_survival_data(data_80, fmt="C1.-", label="8mW")
grid()
legend()
ylim([0, 0.5])
title("N=2 Raman")
xlabel("Raman Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_7-8")

plotx, plote1, plote2, plote3 = eval_model(linspace(1, 9, 1000), fit.param)
uparams = copy(fit.param)
uparams[7:9] = 0
plotx, plotue1, plotue2, plotue3 = eval_model(linspace(1, 9, 1000), uparams)

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

NaCsPlot.maybe_show()
