#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180807_231104.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function gen_selector(na::Bool)
    function selector(logicals)
        @assert size(logicals, 2) == 1
        if logicals[2, 1] == 0
            return Int[1, 0, 0]
        elseif (logicals[1, 1] == 0) == na
            return Int[1, 0, 0]
        end
        return Int[1, 1, logicals[4, 1]]
    end
end

data_cs_a = NaCsData.select_count(params_a, logicals_a, gen_selector(false))
data_nacs_a = NaCsData.select_count(params_a, logicals_a, gen_selector(true))

const spec_a = OrderedDict(
    :raman=>7 .+ linspace(-30, 30, 121),
    :uwave=>7 .+ linspace(-30, 30, 121),
)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)
const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20180807_231104_interaction_shift")

data_cs_raman = split_cs_a[:raman]
data_cs_uwave = split_cs_a[:uwave]
data_nacs_raman = split_nacs_a[:raman]
data_nacs_uwave = split_nacs_a[:uwave]

figure()
NaCsPlot.plot_survival_data(data_cs_raman, fmt="C0.-", label="Cs only")
NaCsPlot.plot_survival_data(data_nacs_raman, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("Raman interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_raman")

figure()
NaCsPlot.plot_survival_data(data_cs_uwave, fmt="C0.-", label="Cs only")
NaCsPlot.plot_survival_data(data_nacs_uwave, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("Microwave interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_uwave")

NaCsPlot.maybe_show()
