#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180808_095637.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function gen_selector(na)
    idx1 = na ? 1 : 2
    idx2 = idx1 + 2
    function selector(logicals)
        @assert size(logicals, 2) == 1
        if logicals[idx1, 1] == 0
            return Int[1, 0, 0]
        end
        return Int[1, 1, logicals[idx2, 1]]
    end
end

data_na_a = NaCsData.select_count(params_a, logicals_a, gen_selector(true))
data_cs_a = NaCsData.select_count(params_a, logicals_a, gen_selector(false))

const spec_a = OrderedDict(
    :p49=>(linspace(-0.8, 0.8, 16), linspace(-0.8, 0.8, 16)),
    :p50=>(linspace(-0.8, 0.8, 16), linspace(-0.8, 0.8, 16)),
    :p51=>(linspace(-0.8, 0.8, 16), linspace(-0.8, 0.8, 16)),
)

const split_na_a = NaCsData.split_data(data_na_a, spec_a)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20180808_095637_merge_heating")

const dMerge = 0.3125

data_na_p49 = NaCsData.map_params((i, v)->v - dMerge, split_na_a[:p49])
data_na_p50 = NaCsData.map_params((i, v)->v, split_na_a[:p50])
data_na_p51 = NaCsData.map_params((i, v)->v + dMerge, split_na_a[:p51])

data_cs_p49 = NaCsData.map_params((i, v)->v - dMerge, split_cs_a[:p49])
data_cs_p50 = NaCsData.map_params((i, v)->v, split_cs_a[:p50])
data_cs_p51 = NaCsData.map_params((i, v)->v + dMerge, split_cs_a[:p51])

figure()
NaCsPlot.plot_survival_data(data_na_p49[1], fmt="C0.-", label="49")
NaCsPlot.plot_survival_data(data_na_p50[1], fmt="C1.-", label="50")
NaCsPlot.plot_survival_data(data_na_p51[1], fmt="C2.-", label="51")
grid()
legend()
ylim([0, 0.6])
title("Na axial (Na=8, Cs=15)")
xlabel("Displacement (um)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_old")

figure()
NaCsPlot.plot_survival_data(data_cs_p49[1], fmt="C0.-", label="49")
NaCsPlot.plot_survival_data(data_cs_p50[1], fmt="C1.-", label="50")
NaCsPlot.plot_survival_data(data_cs_p51[1], fmt="C2.-", label="51")
grid()
legend()
ylim([0, 0.6])
title("Cs axial (Na=8, Cs=15)")
xlabel("Displacement (um)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_old")

figure()
NaCsPlot.plot_survival_data(data_na_p49[2], fmt="C0.-", label="49")
NaCsPlot.plot_survival_data(data_na_p50[2], fmt="C1.-", label="50")
NaCsPlot.plot_survival_data(data_na_p51[2], fmt="C2.-", label="51")
grid()
legend()
ylim([0, 0.6])
title("Na axial (Na=11, Cs=9)")
xlabel("Displacement (um)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_new")

figure()
NaCsPlot.plot_survival_data(data_cs_p49[2], fmt="C0.-", label="49")
NaCsPlot.plot_survival_data(data_cs_p50[2], fmt="C1.-", label="50")
NaCsPlot.plot_survival_data(data_cs_p51[2], fmt="C2.-", label="51")
grid()
legend()
ylim([0, 0.6])
title("Cs axial (Na=11, Cs=9)")
xlabel("Displacement (um)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_new")

NaCsPlot.maybe_show()
