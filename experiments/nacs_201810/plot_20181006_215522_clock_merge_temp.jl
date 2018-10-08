#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181006_215522.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))
data_cs_only_a =
    NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((-1, 2,), (4,)))
data_cs_nacs_a =
    NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2,), (-3, 4)))
data_nacs_nacs_a =
    NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2,), (3, 4)))

const spec_na_a = OrderedDict(
    :cs_shift=>-40:2.0:80,
    :nomerge=>OrderedDict(
        :az=>(-126 .+ (-30:4.0:30), 76 .+ (-30:4.0:30)),
        :rx=>(-479 .+ (-45:6.0:45), 614 .+ (-90:12.0:90)),
        :ry=>(-490 .+ (-45:6.0:45), 631 .+ (-90:12.0:90))
    ),
    :merge=>OrderedDict(
        :az=>(-126 .+ (-30:4.0:30), 76 .+ (-30:4.0:30)),
        :rx=>(-479 .+ (-45:6.0:45), 614 .+ (-90:12.0:90)),
        :ry=>(-490 .+ (-45:6.0:45), 631 .+ (-90:12.0:90))
    ),
)

const spec_cs_a = OrderedDict(
    :cs_shift=>-40:2.0:80,
    :nomerge=>OrderedDict(
        :az=>(-10 .+ (-15:2.0:15), 33 .+ (-15:2.0:15)),
        :rx=>(-135 .+ (-45:6.0:45), 109 .+ (-75:10.0:75)),
        :ry=>(-129 .+ (-45:6.0:45), 112 .+ (-60:8.0:60))
    ),
    :merge=>OrderedDict(
        :az=>(-10 .+ (-15:2.0:15), 33 .+ (-15:2.0:15)),
        :rx=>(-135 .+ (-45:6.0:45), 109 .+ (-75:10.0:75)),
        :ry=>(-129 .+ (-45:6.0:45), 112 .+ (-60:8.0:60))
    ),
)

const split_na_a = NaCsData.split_data(data_na_a, spec_na_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)
const split_cs_only_a = NaCsData.split_data(data_cs_only_a, spec_cs_a)
const split_cs_nacs_a = NaCsData.split_data(data_cs_nacs_a, spec_cs_a)
const split_nacs_nacs_a = NaCsData.split_data(data_nacs_nacs_a, spec_cs_a)

data_cs_only_shift = split_cs_only_a[:cs_shift]
data_cs_nacs_shift = split_cs_nacs_a[:cs_shift]
data_nacs_nacs_shift = split_nacs_nacs_a[:cs_shift]

data_na_nomerge_az = split_na_a[:nomerge][:az]
data_na_nomerge_rx = split_na_a[:nomerge][:rx]
data_na_nomerge_ry = split_na_a[:nomerge][:ry]
data_cs_nomerge_az = split_cs_a[:nomerge][:az]
data_cs_nomerge_rx = split_cs_a[:nomerge][:rx]
data_cs_nomerge_ry = split_cs_a[:nomerge][:ry]

data_na_merge_az = split_na_a[:merge][:az]
data_na_merge_rx = split_na_a[:merge][:rx]
data_na_merge_ry = split_na_a[:merge][:ry]
data_cs_merge_az = split_cs_a[:merge][:az]
data_cs_merge_rx = split_cs_a[:merge][:rx]
data_cs_merge_ry = split_cs_a[:merge][:ry]

const prefix = joinpath(@__DIR__, "imgs", "data_20181006_215522_clock_merge_temp")

figure()
NaCsPlot.plot_survival_data(data_cs_only_shift, fmt="C0.-", label="Cs->Cs")
NaCsPlot.plot_survival_data(data_cs_nacs_shift, fmt="C1.-", label="NaCs->Cs")
NaCsPlot.plot_survival_data(data_nacs_nacs_shift, fmt="C2.-", label="NaCs->NaCs")
grid()
legend()
ylim([0, 1])
title("Cs (4 + 2) -> (3 + 2)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_shift")

figure()
NaCsPlot.plot_survival_data(data_na_nomerge_az[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_na_nomerge_az[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_merge_az[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_na_merge_az[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(data_cs_nomerge_az[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_cs_nomerge_az[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_merge_az[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_cs_merge_az[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(data_na_nomerge_rx[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_na_nomerge_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_merge_rx[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_na_merge_rx[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Na X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
NaCsPlot.plot_survival_data(data_cs_nomerge_rx[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_cs_nomerge_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_merge_rx[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_cs_merge_rx[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Cs X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
NaCsPlot.plot_survival_data(data_na_nomerge_ry[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_na_nomerge_ry[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_merge_ry[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_na_merge_ry[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Na Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
NaCsPlot.plot_survival_data(data_cs_nomerge_ry[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_cs_nomerge_ry[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_merge_ry[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_cs_merge_ry[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Cs Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

NaCsPlot.maybe_show()
