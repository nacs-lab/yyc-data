#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180808_141220.mat")
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

data_cs_a = NaCsData.select_count(params_a, logicals_a, gen_selector(false))

const spec_a = OrderedDict(
    :p15=>[0.1; 0.2; 0.5; 0.8; logspace(0, 1, 51)],
    :p8=>[0.1; 0.2; 0.5; 0.8; logspace(0, 1, 51)],
)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)

data_cs_p8 = split_cs_a[:p8]
data_cs_p15 = split_cs_a[:p15]

const prefix = joinpath(@__DIR__, "imgs", "data_20180808_141220_move_heating")

figure()
NaCsPlot.plot_survival_data(data_cs_p8, fmt="C0.-", label="8")
NaCsPlot.plot_survival_data(data_cs_p15, fmt="C1.-", label="15")
grid()
legend()
ylim([0, 0.6])
title("Cs axial")
xlabel("TMerge (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
