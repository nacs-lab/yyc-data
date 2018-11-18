#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181114_230217.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20181116_113355.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)

data_nacs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2), (3, 4)))
data_nacs_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((1, 2), (3, 4)))

const spec_a = (1.0:(10 * 41), 0.0:0.0)
const split_nacs_a = NaCsData.split_data([data_nacs_a; data_nacs_b], spec_a)

const times = Float64[1, 2, 4, 6, 8, 12, 16, 20, 25, 30]
const freqs = 298.62:0.002:298.70

const time_datas =
    [[split_nacs_a[2];
      NaCsData.map_params((i, v)->times[i],
                          split_nacs_a[1][((i - 1) * 10 + 1):((i - 1) * 10 + 10)])]
     for i in 1:41]
const freq_datas =
    [NaCsData.map_params((i, v)->freqs[i], split_nacs_a[1][i:10:end])
     for i in 1:10]

const prefix = joinpath(@__DIR__, "imgs", "data_20181114_raman_2d")

figure()
for data in time_datas
    NaCsPlot.plot_survival_data(data, fmt=".-")
end
grid()
# legend()
ylim([0, 1])
title("Raman time")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_time")

figure()
for data in freq_datas
    NaCsPlot.plot_survival_data(data, fmt=".-")
end
grid()
# legend()
ylim([0, 1])
title("Raman spectrum")
xlabel("Detuning (\$MHz\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_freq")

NaCsPlot.maybe_show()
