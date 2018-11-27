#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181126_204725.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
data_nacs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2), (3, 4)))

const spec_a = (298.204 .+ (-20:1:20) .* 1e-3,
                298.204 .+ (-20:1:20) .* 1e-3,
                298.204 .+ (-20:1:20) .* 1e-3,
                298.204 .+ (-20:1:20) .* 1e-3)

const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)

data_m10 = [split_nacs_a[1]; split_nacs_a[3]]
data_p10 = [split_nacs_a[2]; split_nacs_a[4]]

const prefix = joinpath(@__DIR__, "imgs", "data_20181126_204725_raman_circ_pol_ang")

figure()
NaCsPlot.plot_survival_data(data_m10, fmt="C0o-", label="-10\${}^\\circ\$")
NaCsPlot.plot_survival_data(data_p10, fmt="C1o-", label="10\${}^\\circ\$")
grid()
legend()
ylim([0.3, 0.8])
title("Raman spectrum")
xlabel("Detuning (\$MHz\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
