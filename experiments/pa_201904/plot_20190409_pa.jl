#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
using MAT

function load_data(names, selector)
    local datas
    for name in names
        param, counts = matopen(joinpath(@__DIR__, "data", "data_$name.mat")) do fd
            sg = read(fd, "ScanGroup")
            f = sg["base"]["params"]["fWavemeter"]
            c = read(fd, "SingleAtomLogical") .!= 0
            return [f for _ in 1:size(c, 3)], c
        end
        data = NaCsData.select_count(param, counts, selector)
        if !@isdefined datas
            datas = data
        else
            datas = [datas; data]
        end
    end
    return datas
end

function load_data3(names, selector)
    local datas
    spec = (1.0:1, 1.0:1, 1.0:1)
    for name in names
        param, counts, freq = matopen(joinpath(@__DIR__, "data", "data_$name.mat")) do fd
            p, c = NaCsData.load_striped_mat(fd)
            sg = read(fd, "ScanGroup")
            f = sg["base"]["params"]["fWavemeter"]
            return p, c, f
        end
        data = NaCsData.map_params((i, v)->freq,
                                   NaCsData.split_data(NaCsData.select_count(param, counts,
                                                                             selector), spec))
        if !@isdefined datas
            datas = data
        else
            datas = ([datas[1]; data[1]], [datas[2]; data[2]], [datas[3]; data[3]])
        end
    end
    return datas
end

const names_files = ["data_20190409_pa_names1.mat",
                     "data_20190409_pa_names2.mat",
                     "data_20190410_pa_names.mat"]
const datas = [load_data(matopen(fd->read(fd, "names"), joinpath(@__DIR__, "data", names_file)),
                         NaCsData.select_single((1, 2,), (3, 4,)))
               for names_file in names_files]

const datas2 = [[NaCsData.map_params((i, v)->((v + 0.05) ÷ 0.1) * 0.1, data);] for data in datas]

const names_files3 = ["data_20190410_pa_names2.mat"]
const datas3 = [load_data3(matopen(fd->read(fd, "names"),
                                   joinpath(@__DIR__, "data", names_file)),
                           NaCsData.select_single((1, 2,), (3, 4,)))
                for names_file in names_files3]

const prefix = joinpath(@__DIR__, "imgs", "data_20190409_pa")

figure()
NaCsPlot.plot_survival_data(datas[1], fmt="C0.-")
NaCsPlot.plot_survival_data(datas[2], fmt="C0.-")
NaCsPlot.plot_survival_data(datas[3], fmt="C0.-")
grid()
title("PA spectrum")
xlabel("288XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

figure()
NaCsPlot.plot_survival_data(datas2[1], fmt="C0.-")
NaCsPlot.plot_survival_data(datas2[2], fmt="C0.-")
NaCsPlot.plot_survival_data(datas2[3], fmt="C0.-")
grid()
title("PA spectrum (averaged)")
xlabel("288XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_avg")

figure()
NaCsPlot.plot_survival_data(datas2[1], fmt="C0.-", label="0.8ms 5.8mW")
NaCsPlot.plot_survival_data(datas3[1][1], fmt="C1.-", label="1.0ms 0.53mW")
NaCsPlot.plot_survival_data(datas3[1][2], fmt="C2.-", label="0.8ms 1.0mW")
NaCsPlot.plot_survival_data(datas3[1][3], fmt="C3.-", label="0.3ms 5.8mW")
ylim([0.2, 0.9])
xlim([692.8, 695.2])
legend(fontsize="small", ncol=2, labelspacing=0.2, borderpad=0.4,
       handletextpad=0.2, columnspacing=0.4, borderaxespad=0.2)
grid()
title("PA spectrum @ ms, mW")
xlabel("288XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)3")

NaCsPlot.maybe_show()
