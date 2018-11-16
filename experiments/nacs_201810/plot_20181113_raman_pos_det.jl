#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181112_153312.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20181113_232614.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const iname_c = joinpath(@__DIR__, "data", "data_20181114_074540.mat")
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)
const iname_d = joinpath(@__DIR__, "data", "data_20181113_081215.mat")
const params_d, logicals_d = NaCsData.load_striped_mat(iname_d)
const iname_e = joinpath(@__DIR__, "data", "data_20181113_095954.mat")
const params_e, logicals_e = NaCsData.load_striped_mat(iname_e)

data_nacs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2), (3, 4)))
data_nacs_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((1, 2), (3, 4)))
data_nacs_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((1, 2), (3, 4)))
data_nacs_d = NaCsData.select_count(params_d, logicals_d, NaCsData.select_single((1, 2), (3, 4)))
data_nacs_e = NaCsData.select_count(params_e, logicals_e, NaCsData.select_single((1, 2), (3, 4)))

const spec_a = linspace(0, 8, 41)
const spec_b = (293:0.05:294.5, linspace(0, 8, 41))
const spec_c = linspace(0, 8, 41)
const spec_d = (301.2:0.05:302.5, linspace(0, 8, 41), linspace(0, 8, 41))
const spec_e = (linspace(0, 8, 41), linspace(0, 8, 41))

const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)
const split_nacs_b = NaCsData.split_data(data_nacs_b, spec_b)
const split_nacs_c = NaCsData.split_data(data_nacs_c, spec_c)
const split_nacs_d = NaCsData.split_data(data_nacs_d, spec_d)
const split_nacs_e = NaCsData.split_data(data_nacs_e, spec_e)

data_time_m20 = [split_nacs_a; split_nacs_d[2]; split_nacs_e[1]]
data_freq_m20 = split_nacs_d[1]
data_time_p20 = [split_nacs_b[2]; split_nacs_c]
data_freq_p20 = split_nacs_b[1]

const prefix = joinpath(@__DIR__, "imgs", "data_20181113_raman_pos_det")

figure()
NaCsPlot.plot_survival_data(data_time_m20, fmt="C0.-", label="\$-20GHz\$")
NaCsPlot.plot_survival_data(data_time_p20, fmt="C1.-", label="\$+20GHz\$")
grid()
legend()
ylim([0, 1])
title("Raman time")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_time")

figure()
NaCsPlot.plot_survival_data(data_freq_m20, fmt="C0.-", label="\$-20GHz\$")
NaCsPlot.plot_survival_data(data_freq_p20, fmt="C1.-", label="\$+20GHz\$")
grid()
legend()
ylim([0, 1])
title("Raman spectrum")
xlabel("Detuning (\$MHz\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_freq")

NaCsPlot.maybe_show()
