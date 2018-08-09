#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180808_174308.mat")
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
data_na_a = NaCsData.select_count(params_a, logicals_a, gen_selector(true))

const spec_a = OrderedDict(
    :p15=>linspace(1, 30, 30),
    :p9=>linspace(1, 30, 30),
)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)
const split_na_a = NaCsData.split_data(data_na_a, spec_a)

data_cs_p9 = split_cs_a[:p9]
data_cs_p15 = split_cs_a[:p15]

data_na_p9 = split_na_a[:p9]
data_na_p15 = split_na_a[:p15]

const prefix = joinpath(@__DIR__, "imgs", "data_20180808_174308_napremerge")

figure()
NaCsPlot.plot_survival_data(data_cs_p15, fmt="C0.-", label="Cs")
NaCsPlot.plot_survival_data(data_na_p15, fmt="C1.-", label="Na")
grid()
legend()
ylim([0, 0.6])
title("Axial (Cs=15mW)")
xlabel("NaPreMerge (mW)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs15")

figure()
NaCsPlot.plot_survival_data(data_cs_p9, fmt="C0.-", label="Cs")
NaCsPlot.plot_survival_data(data_na_p9, fmt="C1.-", label="Na")
grid()
legend()
ylim([0, 0.6])
title("Axial (Cs=9mW)")
xlabel("NaPreMerge (mW)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs9")

NaCsPlot.maybe_show()
