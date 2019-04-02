#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190321_234645.mat",
                "data_20190322_093053.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int), typemax(Int)]
const specs_na = [(OrderedDict(
    :az=>(-130 .+ (-30:4.0:30), 67 .+ (-30:4.0:30)),
    :rx=>(-480 .+ (-45:6.0:45), 549 .+ (-90:12.0:90)),
    :ry=>(-481 .+ (-45:6.0:45), 553 .+ (-90:12.0:90))
), OrderedDict(
    :az=>(-130 .+ (-30:4.0:30), 67 .+ (-30:4.0:30)),
    :rx=>(-480 .+ (-45:6.0:45), 549 .+ (-90:12.0:90)),
    :ry=>(-481 .+ (-45:6.0:45), 553 .+ (-90:12.0:90))
), -80:4:60, -80:4:60, -80:4:60),
                  -80:4:60]
const specs_cs = [(OrderedDict(
    :az=>(-28 .+ (-15:2.0:15), 14 .+ (-15:2.0:15)),
    :rx=>(-155 .+ (-45:6.0:45), 85 .+ (-75:10.0:75)),
    :ry=>(-160 .+ (-45:6.0:45), 80 .+ (-60:8.0:60))
), OrderedDict(
    :az=>(-28 .+ (-15:2.0:15), 14 .+ (-15:2.0:15)),
    :rx=>(-155 .+ (-45:6.0:45), 85 .+ (-75:10.0:75)),
    :ry=>(-160 .+ (-45:6.0:45), 80 .+ (-60:8.0:60))
), -80:4:60, -80:4:60, -80:4:60),
                  -80:4:60]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]
const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)
const datas_nacs_na = select_datas(datas, NaCsData.select_single((1, 2), (3,)), maxcnts, specs_na)
const datas_nacs_cs = select_datas(datas, NaCsData.select_single((1, 2), (4,)), maxcnts, specs_cs)
const datas_cs_cs = select_datas(datas, NaCsData.select_single((-1, 2), (4,)), maxcnts, specs_cs)

data_na_az = datas_na[1][1][:az]
data_na_rx = datas_na[1][1][:rx]
data_na_ry = datas_na[1][1][:ry]
data_cs_az = datas_cs[1][1][:az]
data_cs_rx = datas_cs[1][1][:rx]
data_cs_ry = datas_cs[1][1][:ry]
data_na_az_merge = datas_na[1][2][:az]
data_na_rx_merge = datas_na[1][2][:rx]
data_na_ry_merge = datas_na[1][2][:ry]
data_cs_az_merge = datas_cs[1][2][:az]
data_cs_rx_merge = datas_cs[1][2][:rx]
data_cs_ry_merge = datas_cs[1][2][:ry]
data_na_shift = [datas_nacs_na[1][3]; datas_nacs_na[1][4]; datas_nacs_na[1][5]; datas_nacs_na[2]]
data_cs_shift = [datas_nacs_cs[1][3]; datas_nacs_cs[1][4]; datas_nacs_cs[1][5]; datas_nacs_cs[2]]
data_cs_cs_shift = [datas_cs_cs[1][3]; datas_cs_cs[1][4]; datas_cs_cs[1][5]; datas_cs_cs[2]]

const prefix = joinpath(@__DIR__, "imgs", "data_20190321_temp")
const prefix2 = joinpath(@__DIR__, "imgs", "data_20190321_clock_shift")

figure()
NaCsPlot.plot_survival_data(data_na_az[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_na_az[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_az_merge[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_na_az_merge[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(data_cs_az[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_cs_az[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_az_merge[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_cs_az_merge[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(data_na_rx[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_na_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_rx_merge[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_na_rx_merge[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Na X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
NaCsPlot.plot_survival_data(data_cs_rx[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_cs_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_rx_merge[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_cs_rx_merge[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Cs X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
NaCsPlot.plot_survival_data(data_na_ry[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_na_ry[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_ry_merge[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_na_ry_merge[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Na Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
NaCsPlot.plot_survival_data(data_cs_ry[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_cs_ry[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_ry_merge[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_cs_ry_merge[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Cs Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

figure()
NaCsPlot.plot_survival_data(data_na_shift, fmt="C0.-", label="NaX / NaCs")
NaCsPlot.plot_survival_data(data_cs_shift, fmt="C1.-", label="CsX / NaCs")
NaCsPlot.plot_survival_data(data_cs_cs_shift, fmt="C2.-", label="Cs / Cs")
legend(fontsize="small")
grid()
ylim([0, 1])
title("Interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix2)")

NaCsPlot.maybe_show()
