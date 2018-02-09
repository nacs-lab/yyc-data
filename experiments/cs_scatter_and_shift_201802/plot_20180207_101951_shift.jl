#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Atomic
import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180207_101951.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180207_180959.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)

function gen_selector(na)
    function selector(logicals)
        @assert size(logicals, 2) == 1
        if logicals[2, 1] == 0 || logicals[1, 1] != na
            return Int[1, 0, 0]
        end
        return Int[1, 1, logicals[4, 1]]
    end
end

data_single_a = NaCsData.select_count(params_a, logicals_a, gen_selector(false))
data_both_a = NaCsData.select_count(params_a, logicals_a, gen_selector(true))
data_single_b = NaCsData.select_count(params_b, logicals_b, gen_selector(false))
data_both_b = NaCsData.select_count(params_b, logicals_b, gen_selector(true))

const spec_a = OrderedDict(
    :set1=>linspace(-50, 130, 61),
    :set2=>linspace(-50, 130, 61),
)

const spec_b = OrderedDict(
    :set1=>linspace(-50, 130, 31),
    :set2=>linspace(-50, 130, 31),
)

const split_single_a = NaCsData.split_data(data_single_a, spec_a)
const split_both_a = NaCsData.split_data(data_both_a, spec_a)
const split_single_b = NaCsData.split_data(data_single_b, spec_b)
const split_both_b = NaCsData.split_data(data_both_b, spec_b)

const data_single = [split_single_a[:set1]; split_single_a[:set2];
                     split_single_b[:set1]; split_single_b[:set2]]
const data_both = [split_both_a[:set1]; split_both_a[:set2];
                   split_both_b[:set1]; split_both_b[:set2]]

const prefix = joinpath(@__DIR__, "imgs", "data_20180207_101951_shift")

figure()
NaCsPlot.plot_survival_data(data_single, fmt="C0o-", label="Cs only")
NaCsPlot.plot_survival_data(data_both, fmt="C1o-", label="Cs + Na")
grid()
ylim([0, 1])
legend()
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
