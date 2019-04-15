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

const names_files = ["data_20190414_pa_names2.mat"]
const datas = [load_data(matopen(fd->read(fd, "names"), joinpath(@__DIR__, "data", names_file)),
                         NaCsData.select_single((1, 2,), (3, 4,)))
               for names_file in names_files]
const datas_na = [load_data(matopen(fd->read(fd, "names"),
                                    joinpath(@__DIR__, "data", names_file)),
                            NaCsData.select_single((1, 2,), (3,)))
                  for names_file in names_files]
const datas_cs = [load_data(matopen(fd->read(fd, "names"),
                                    joinpath(@__DIR__, "data", names_file)),
                            NaCsData.select_single((1, 2,), (4,)))
                  for names_file in names_files]
const datas_cs_cs = [load_data(matopen(fd->read(fd, "names"),
                                       joinpath(@__DIR__, "data", names_file)),
                               NaCsData.select_single((-1, 2,), (-3, 4,)))
                     for names_file in names_files]

const prefix = joinpath(@__DIR__, "imgs", "data_20190414_pa")

figure()
NaCsPlot.plot_survival_data(datas[1], fmt="C0.-", label="Na+Cs")
NaCsPlot.plot_survival_data(datas_na[1], fmt="C1.-", label="Na")
NaCsPlot.plot_survival_data(datas_cs[1], fmt="C2.-", label="Cs")
legend(fontsize="small", ncol=3)
ylim([0.35, 0.9])
grid()
title("PA spectrum")
xlabel("288XXX GHz")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

figure()
NaCsPlot.plot_survival_data(datas_cs[1], fmt="C0.-", label="Cs/Na+Cs")
NaCsPlot.plot_survival_data(datas_cs_cs[1], fmt="C1.-", label="Cs/Cs")
legend(fontsize="small", ncol=2)
ylim([0.35, 0.9])
grid()
title("Cs survival")
xlabel("288XXX GHz")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
