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

const names_files = ["data_20190603_pa32_names.mat"]
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
const datas_na_na = [load_data(matopen(fd->read(fd, "names"),
                                       joinpath(@__DIR__, "data", names_file)),
                               NaCsData.select_single((1, -2,), (3, -4,)))
                     for names_file in names_files]

const prefix = joinpath(@__DIR__, "imgs", "data_20190603_pa32")

figure()
NaCsPlot.plot_survival_data(datas[1], fmt="C0.-")
grid()
title("PA spectrum")
xlabel("307XXX GHz")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
