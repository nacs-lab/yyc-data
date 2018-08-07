#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180806_104858.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180806_153059.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const iname_c = joinpath(@__DIR__, "data", "data_20180806_173202.mat")
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)
const iname_d = joinpath(@__DIR__, "data", "data_20180806_203632.mat")
const params_d, logicals_d = NaCsData.load_striped_mat(iname_d)

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
data_na_b = NaCsData.select_count(params_b, logicals_b, gen_selector(true))
data_na_c = NaCsData.select_count(params_c, logicals_c, gen_selector(true))
data_na_d = NaCsData.select_count(params_d, logicals_d, gen_selector(true))

data_cs_a = NaCsData.select_count(params_a, logicals_a, gen_selector(false))
data_cs_b = NaCsData.select_count(params_b, logicals_b, gen_selector(false))
data_cs_c = NaCsData.select_count(params_c, logicals_c, gen_selector(false))
data_cs_d = NaCsData.select_count(params_d, logicals_d, gen_selector(false))

# Na hot +1, -1
const spec_na_a = OrderedDict(
    :az=>(-109 .+ 30 .* linspace(-1, 1, 16), 64 .+ 30 .* linspace(-1, 1, 16)),
    :rx=>(-421 .+ 45 .* linspace(-1, 1, 16), 579 .+ 90 .* linspace(-1, 1, 16)),
    :ry=>(-430 .+ 45 .* linspace(-1, 1, 16), 588 .+ 90 .* linspace(-1, 1, 16)),
)
# Cs bad cooling +1, -1
const spec_cs_a = OrderedDict(
    :az=>(-21 .+ 15 .* linspace(-1, 1, 16), 20 .+ 15 .* linspace(-1, 1, 16)),
    :rx=>(-153 .+ 45 .* linspace(-1, 1, 16), 100 .+ 75 .* linspace(-1, 1, 16)),
    :ry=>(-149 .+ 45 .* linspace(-1, 1, 16), 110 .+ 60 .* linspace(-1, 1, 16)),
)

# Na hot ax -2, -3, rad -2
const spec_na_b = OrderedDict(
    :az=>(136 .+ 30 .* linspace(-1, 1, 16), 224 .+ 30 .* linspace(-1, 1, 16)),
    :rx=>(1000 .+ 90 .* linspace(-1, 1, 16)),
    :ry=>(1068 .+ 90 .* linspace(-1, 1, 16)),
)
# Cs bad cooling ax +1, -1 higher power, rad -1
const spec_cs_b = OrderedDict(
    :az=>(-21 .+ 15 .* linspace(-1, 1, 16), 20 .+ 15 .* linspace(-1, 1, 16)), # higher power
    :rx=>(100 .+ 75 .* linspace(-1, 1, 16)),
    :ry=>(110 .+ 60 .* linspace(-1, 1, 16)),
)

# Na hot ax -4, -5
const spec_na_c = (301 .+ 30 .* linspace(-1, 1, 16), 384 .+ 30 .* linspace(-1, 1, 16))
# Cs normal cooling ax +1, -1
const spec_cs_c = (-21 .+ 15 .* linspace(-1, 1, 16), 20 .+ 15 .* linspace(-1, 1, 16))

# Na hot rabi -1
const spec_na_d = OrderedDict(
    :az=>linspace(0, 180, 21),
    :rx=>linspace(0, 80, 21),
    :ry=>linspace(0, 80, 21),
)
# Cs normal rabi +1
const spec_cs_d = OrderedDict(
    :az=>linspace(0, 800, 21),
    :rx=>linspace(0, 240, 21),
    :ry=>linspace(0, 240, 21),
)

const split_na_a = NaCsData.split_data(data_na_a, spec_na_a)
const split_na_b = NaCsData.split_data(data_na_b, spec_na_b)
const split_na_c = NaCsData.split_data(NaCsData.map_params((i, f)->i, data_na_c), spec_na_c)
const split_na_d = NaCsData.split_data(data_na_d, spec_na_d)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)
const split_cs_b = NaCsData.split_data(data_cs_b, spec_cs_b)
const split_cs_c = NaCsData.split_data(NaCsData.map_params((i, f)->i, data_cs_c), spec_cs_c)
const split_cs_d = NaCsData.split_data(data_cs_d, spec_cs_d)

const prefix = joinpath(@__DIR__, "imgs", "data_20180806_hot")

figure()
NaCsPlot.plot_survival_data(split_na_a[:az][1], fmt="C0.-") # -119
NaCsPlot.plot_survival_data(split_na_a[:az][2], fmt="C0.-") # 50
NaCsPlot.plot_survival_data(split_na_b[:az][1], fmt="C0.-") # 136
NaCsPlot.plot_survival_data(split_na_b[:az][2], fmt="C0.-") # 222
NaCsPlot.plot_survival_data(split_na_c[1], fmt="C0.-") # 303
NaCsPlot.plot_survival_data(split_na_c[2], fmt="C0.-") # 388
grid()
ylim([0, 0.6])
title("Na Hot Axis Z (axial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(split_na_a[:rx][1], fmt="C0.-") # -430
NaCsPlot.plot_survival_data(split_na_a[:rx][2], fmt="C0.-") # 560
NaCsPlot.plot_survival_data(split_na_b[:rx], fmt="C0.-") # 1035
grid()
ylim([0, 0.65])
title("Na Hot Axis X (radial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
NaCsPlot.plot_survival_data(split_na_a[:ry][1], fmt="C0.-") # -442
NaCsPlot.plot_survival_data(split_na_a[:ry][2], fmt="C0.-") # 568
NaCsPlot.plot_survival_data(split_na_b[:ry], fmt="C0.-") # 1040
grid()
ylim([0, 0.7])
title("Na Hot Axis Y (radial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
NaCsPlot.plot_survival_data(split_cs_a[:az][1], fmt="C0.-", label="Narrow")
NaCsPlot.plot_survival_data(split_cs_a[:az][2], fmt="C0.-")
NaCsPlot.plot_survival_data(split_cs_b[:az][1], fmt="C1.-", label="Wide") # -20
NaCsPlot.plot_survival_data(split_cs_b[:az][2], fmt="C1.-") # 21
NaCsPlot.plot_survival_data(split_cs_c[1], fmt="C2.-", label="Wide 2")
NaCsPlot.plot_survival_data(split_cs_c[2], fmt="C2.-")
grid()
legend()
ylim([0, 0.8])
title("Cs Axis Z (axial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(split_cs_a[:rx][1], fmt="C0.-") # -153
NaCsPlot.plot_survival_data([split_cs_a[:rx][2]; split_cs_b[:rx]], fmt="C0.-") # 95
grid()
ylim([0, 0.85])
title("Cs Hot Axis X (radial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
NaCsPlot.plot_survival_data(split_cs_a[:ry][1], fmt="C0.-") # -148
NaCsPlot.plot_survival_data([split_cs_a[:ry][2]; split_cs_b[:ry]], fmt="C0.-") # 119
grid()
ylim([0, 0.85])
title("Cs Hot Axis Y (radial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

figure()
NaCsPlot.plot_survival_data(split_na_d[:az], fmt="C0.-", label="Axial Z")
NaCsPlot.plot_survival_data(split_na_d[:rx], fmt="C1.-", label="Radial X")
NaCsPlot.plot_survival_data(split_na_d[:ry], fmt="C2.-", label="Radial Y")
grid()
legend()
ylim([0, 0.65])
title("Na Hot -1 Rabi flopping (axial)")
xlabel("Time (us)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rabi")

figure()
NaCsPlot.plot_survival_data(split_cs_d[:az], fmt="C0.-", label="Axial Z")
NaCsPlot.plot_survival_data(split_cs_d[:rx], fmt="C1.-", label="Radial X")
NaCsPlot.plot_survival_data(split_cs_d[:ry], fmt="C2.-", label="Radial Y")
grid()
legend()
ylim([0, 0.9])
title("Cs +1 Rabi flopping (axial)")
xlabel("Time (us)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rabi")

NaCsPlot.maybe_show()
