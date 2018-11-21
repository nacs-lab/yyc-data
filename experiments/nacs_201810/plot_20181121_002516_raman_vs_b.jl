#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181121_002516.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_nacs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2), (3, 4)))

const spec_a = (298.138 .+ (-5:1:5) .* 1e-3,
                298.138 .+ (-5:1:5) .* 1e-3,
                298.138 .+ (-5:1:5) .* 1e-3)

const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)

data_norm = split_nacs_a[1]
data_b60 = split_nacs_a[2]
data_b60_45 = split_nacs_a[3]

const prefix = joinpath(@__DIR__, "imgs", "data_20181121_002516_raman_vs_b")

figure()
NaCsPlot.plot_survival_data(data_norm, fmt="C0.-", label="100%")
NaCsPlot.plot_survival_data(data_b60, fmt="C1.-", label="60%")
NaCsPlot.plot_survival_data(data_b60_45, fmt="C2.-", label="60%, 45\${}^\\circ\$")
grid()
legend()
ylim([0.4, 0.7])
title("Raman spectrum vs B field")
xlabel("Detuning (\$MHz\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
