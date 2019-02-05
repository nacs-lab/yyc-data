#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190203_094225.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [(OrderedDict(
    :az=>(65 .+ (-30:4.0:30)),
    :rx=>(555 .+ (-90:12.0:90)),
    :ry=>(559 .+ (-90:12.0:90))
), OrderedDict(
    :az=>(65 .+ (-30:4.0:30)),
    :rx=>(555 .+ (-90:12.0:90)),
    :ry=>(559 .+ (-90:12.0:90))
))]
const specs_cs = [(OrderedDict(
    :az=>(16 .+ (-15:2.0:15)),
    :rx=>(89 .+ (-75:10.0:75)),
    :ry=>(112 .+ (-60:8.0:60))
), OrderedDict(
    :az=>(16 .+ (-15:2.0:15)),
    :rx=>(89 .+ (-75:10.0:75)),
    :ry=>(112 .+ (-60:8.0:60))
))]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]
const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_na_az = (datas_na[1][1][:az], datas_na[1][2][:az])
data_na_rx = (datas_na[1][1][:rx], datas_na[1][2][:rx])
data_na_ry = (datas_na[1][1][:ry], datas_na[1][2][:ry])
data_cs_az = [datas_cs[1][1][:az]; datas_cs[1][2][:az]]
data_cs_rx = [datas_cs[1][1][:rx]; datas_cs[1][2][:rx]]
data_cs_ry = [datas_cs[1][1][:ry]; datas_cs[1][2][:ry]]

const prefix = joinpath(@__DIR__, "imgs", "data_20190203_094225_temp")

function freq_arrow(x, y, len, color)
    xl = xlim()
    arrow(x, y, 0, -len, color=color,
          width=(xl[2] - xl[1]) / 80, length_includes_head=true, head_length=len / 3)
end

figure()
NaCsPlot.plot_survival_data(data_na_az[1], fmt="C0.-", label="X555, Y559")
NaCsPlot.plot_survival_data(data_na_az[2], fmt="C1.-", label="X585, Y577")
legend()
grid()
ylim([0, 0.2])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(data_na_rx[1], fmt="C0.-", label="X555, Y559")
NaCsPlot.plot_survival_data(data_na_rx[2], fmt="C1.-", label="X585, Y577")
legend()
grid()
freq_arrow(555, 0.115, 0.02, "C0")
freq_arrow(585, 0.13, 0.02, "C1")
ylim([0, 0.2])
title("Na X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
NaCsPlot.plot_survival_data(data_na_ry[1], fmt="C0.-", label="X555, Y559")
NaCsPlot.plot_survival_data(data_na_ry[2], fmt="C1.-", label="X585, Y577")
legend()
grid()
freq_arrow(559, 0.08, 0.02, "C0")
freq_arrow(577, 0.08, 0.02, "C1")
ylim([0, 0.2])
title("Na Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
NaCsPlot.plot_survival_data(data_cs_az, fmt="C0.-")
grid()
ylim([0, 0.2])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(data_cs_rx, fmt="C0.-")
grid()
ylim([0, 0.2])
title("Cs X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
NaCsPlot.plot_survival_data(data_cs_ry, fmt="C0.-")
grid()
ylim([0, 0.2])
title("Cs Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

NaCsPlot.maybe_show()
