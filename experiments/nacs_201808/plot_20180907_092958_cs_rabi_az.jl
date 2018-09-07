#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180907_092958.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_cs_a = OrderedDict(
    :nm1=>0:40.0:800,
    :n0=>0:40.0:800,
)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20180907_092958_cs_rabi_az")

figure()
NaCsPlot.plot_survival_data(split_cs_a[:nm1], fmt="C0.-", label="+1 order")
NaCsPlot.plot_survival_data(split_cs_a[:n0], fmt="C1.-", label="0 order")
grid()
legend()
ylim([0, 0.9])
title("Cs Axial Rabi flopping")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
