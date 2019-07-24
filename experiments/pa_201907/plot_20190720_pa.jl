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

const names_files = ["data_20190720_pa_names.mat",
                     "data_20190721_pa_names.mat"]
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

const prefix = joinpath(@__DIR__, "imgs", "data_20190720_pa")

figure()
NaCsPlot.plot_survival_data([datas_cs[1]; datas_cs[2]], fmt="C2.-", label="Cs/NaCs")
NaCsPlot.plot_survival_data([datas_na[1]; datas_na[2]], fmt="C1.-", label="Na/NaCs")
NaCsPlot.plot_survival_data([datas[1]; datas[2]], fmt="C0.-", label="NaCs/NaCs")
grid()
ylim([0.35, 0.95])
legend(fontsize="small")
title("PA spectrum (two-body)")
xlabel("307XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_2b")

figure()
NaCsPlot.plot_survival_data([datas_cscs[1]; datas_cscs[2]], fmt="C2.-", label="Cs")
NaCsPlot.plot_survival_data([datas_nana[1]; datas_nana[2]], fmt="C1.-", label="Na")
grid()
ylim([0.35, 0.95])
legend(fontsize="small")
title("PA spectrum (one-body)")
xlabel("307XXX GHz")
ylabel("One-body survival")
NaCsPlot.maybe_save("$(prefix)_1b")

NaCsPlot.maybe_show()
