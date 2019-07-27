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
            pl = read(fd, "ParamList")
            sg = read(fd, "ScanGroup")
            f = sg["base"]["params"]["fWavemeter"]
            amps = sg["base"]["vars"]["params"]["Merge"]["PA"]["Amp"]
            c = read(fd, "SingleAtomLogical") .!= 0
            return [(f, amps[Int(pl[i])]) for i in 1:size(c, 3)], c
        end
        data = NaCsData.select_count(param, counts, selector)
        if !@isdefined datas
            datas = data
        else
            datas = [datas; data]
        end
    end
    datas = datas[[1:2:size(datas, 1) 2:2:size(datas, 1)]]
    datas = NaCsData.map_params((i1, i2, v)->v[1], datas)
    return datas
end

const names_files = ["names_20190725_135756.mat"]
const datas = [load_data(matopen(fd->read(fd, "names"), joinpath(@__DIR__, "data", names_file)),
                         NaCsData.select_single((1, 2,), (3, 4,)))
               for names_file in names_files]
const datas_na = [load_data(matopen(fd->read(fd, "names"), joinpath(@__DIR__, "data", names_file)),
                            NaCsData.select_single((1, 2,), (3,)))
                  for names_file in names_files]
const datas_cs = [load_data(matopen(fd->read(fd, "names"), joinpath(@__DIR__, "data", names_file)),
                            NaCsData.select_single((1, 2,), (4,)))
                  for names_file in names_files]
const datas_nana = [load_data(matopen(fd->read(fd, "names"), joinpath(@__DIR__, "data", names_file)),
                              NaCsData.select_single((1, -2,), (3,)))
                    for names_file in names_files]
const datas_cscs = [load_data(matopen(fd->read(fd, "names"), joinpath(@__DIR__, "data", names_file)),
                              NaCsData.select_single((-1, 2,), (4,)))
                    for names_file in names_files]

const prefix = joinpath(@__DIR__, "imgs", "data_20190725_135756_pa")

figure()
NaCsPlot.plot_survival_data(datas[1][:, 1], fmt="C0.-", label="0.3")
NaCsPlot.plot_survival_data(datas[1][:, 2], fmt="C1.-", label="0.75")
grid()
ylim([0.55, 0.9])
legend(fontsize="small", ncol=2)
title("PA spectrum (NaCs/NaCs)")
xlabel("307XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_2b")

figure()
NaCsPlot.plot_survival_data(datas_cs[1][:, 1], fmt="C2.-", label="Cs/NaCs")
NaCsPlot.plot_survival_data(datas_na[1][:, 1], fmt="C1.-", label="Na/NaCs")
NaCsPlot.plot_survival_data(datas[1][:, 1], fmt="C0.-", label="NaCs/NaCs")
grid()
ylim([0.55, 0.9])
legend(fontsize="small", ncol=2)
title("PA spectrum (two-body 0.30)")
xlabel("307XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_2b_30")

figure()
NaCsPlot.plot_survival_data(datas_cs[1][:, 2], fmt="C2.-", label="Cs/NaCs")
NaCsPlot.plot_survival_data(datas_na[1][:, 2], fmt="C1.-", label="Na/NaCs")
NaCsPlot.plot_survival_data(datas[1][:, 2], fmt="C0.-", label="NaCs/NaCs")
grid()
ylim([0.55, 0.9])
legend(fontsize="small", ncol=2)
title("PA spectrum (two-body 0.75)")
xlabel("307XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_2b_75")

figure()
NaCsPlot.plot_survival_data(datas_cscs[1][:, 1], fmt="C2.-", label="Cs 0.3")
NaCsPlot.plot_survival_data(datas_nana[1][:, 1], fmt="C1.-", label="Na 0.3")
NaCsPlot.plot_survival_data(datas_cscs[1][:, 2], fmt="C2.--", label="Cs 0.75")
NaCsPlot.plot_survival_data(datas_nana[1][:, 2], fmt="C1.--", label="Na 0.75")
grid()
ylim([0.55, 0.9])
legend(fontsize="small", ncol=2)
title("PA spectrum (one-body)")
xlabel("307XXX GHz")
ylabel("One-body survival")
NaCsPlot.maybe_save("$(prefix)_1b")

NaCsPlot.maybe_show()
