#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180815_000658.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180815_091610.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const iname_c = joinpath(@__DIR__, "data", "data_20180815_231138.mat")
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)

data_cs_a = NaCsData.select_count(params_a, logicals_a,
                                  NaCsData.select_single((-1, 2), (4,)))
data_nacs_a = NaCsData.select_count(params_a, logicals_a,
                                    NaCsData.select_single((1, 2), (4,)))
data_cs_b = NaCsData.select_count(params_b, logicals_b,
                                  NaCsData.select_single((-1, 2), (4,)))
data_nacs_b = NaCsData.select_count(params_b, logicals_b,
                                    NaCsData.select_single((1, 2), (4,)))
data_cs_c = NaCsData.select_count(params_c, logicals_c,
                                  NaCsData.select_single((-1, 2), (4,)))
data_nacs_cs_c = NaCsData.select_count(params_c, logicals_c,
                                       NaCsData.select_single((1, 2), (4,)))
data_na_c = NaCsData.select_count(params_c, logicals_c,
                                  NaCsData.select_single((1, -2), (3,)))
data_nacs_na_c = NaCsData.select_count(params_c, logicals_c,
                                       NaCsData.select_single((1, 2), (3,)))

const spec_a = OrderedDict(
    :n0=>7 .+ linspace(-40, 40, 81),
    :n1=>7 .+ linspace(-40, 40, 81),
)
const spec_b = OrderedDict(
    :n0=>7 .+ linspace(-60.0, 10.0, 81),
    :n1=>7 .+ linspace(-60.0, 10.0, 81),
)
const spec_c = OrderedDict(
    :cs=>linspace(-80.0, 80.0, 161),
    :na=>linspace(-80.0, 80.0, 161),
)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)
const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)
const split_cs_b = NaCsData.split_data(data_cs_b, spec_b)
const split_nacs_b = NaCsData.split_data(data_nacs_b, spec_b)
const split_cs_c = NaCsData.split_data(data_cs_c, spec_c)
const split_nacs_cs_c = NaCsData.split_data(data_nacs_cs_c, spec_c)
const split_na_c = NaCsData.split_data(data_na_c, spec_c)
const split_nacs_na_c = NaCsData.split_data(data_nacs_na_c, spec_c)

const prefix = joinpath(@__DIR__, "imgs", "data_20180815_interaction_shift4")

data_cs_n0 = [split_cs_a[:n0]; split_cs_b[:n0]; split_cs_c[:cs]]
data_cs_n1 = [split_cs_a[:n1]; split_cs_b[:n1]]
data_nacs_cs_n0 = [split_nacs_a[:n0]; split_nacs_b[:n0]; split_nacs_cs_c[:cs]]
data_nacs_cs_n1 = [split_nacs_a[:n1]; split_nacs_b[:n1]]

data_na_n0 = split_na_c[:na]
data_nacs_na_n0 = split_nacs_na_c[:na]

figure()
NaCsPlot.plot_survival_data(data_na_n0, fmt="C0.-", label="Na only")
NaCsPlot.plot_survival_data(data_nacs_na_n0, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("(3+2 -> 3+1) n=0 interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_n0")

figure()
NaCsPlot.plot_survival_data(data_cs_n0, fmt="C0.-", label="Cs only")
NaCsPlot.plot_survival_data(data_nacs_cs_n0, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("(4+2 -> 3+2) n=0 interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_n0")

figure()
NaCsPlot.plot_survival_data(data_cs_n1, fmt="C0.-", label="Cs only")
NaCsPlot.plot_survival_data(data_nacs_cs_n1, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("(4+2 -> 3+2) n=1 interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_n1")

figure()
NaCsPlot.plot_survival_data(data_cs_n0, fmt="C0.-", label="n=0")
NaCsPlot.plot_survival_data(data_cs_n1, fmt="C1.-", label="n=1")
grid()
legend()
ylim([0, 0.8])
title("(4+2 -> 3+2) Cs only")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

figure()
NaCsPlot.plot_survival_data(data_nacs_cs_n0, fmt="C0.-", label="n=0")
NaCsPlot.plot_survival_data(data_nacs_cs_n1, fmt="C1.-", label="n=1")
grid()
legend()
ylim([0, 0.8])
title("(4+2 -> 3+2) Na + Cs")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs_cs")

NaCsPlot.maybe_show()
