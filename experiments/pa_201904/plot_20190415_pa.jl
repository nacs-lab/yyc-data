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

const names_files = ["data_20190415_pa_names.mat"]
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

const prefix = joinpath(@__DIR__, "imgs", "data_20190415_pa")

figure()
NaCsPlot.plot_survival_data(datas[1][1], fmt="C0.-", label="\$\\pi\$")
NaCsPlot.plot_survival_data(datas[1][2], fmt="C1.-", label="\$\\sigma\$")
NaCsPlot.plot_survival_data(datas[1][3], fmt="C2.-", label="\$\\pi\$, 10mW")
legend(fontsize="small", ncol=3)
ylim([0.3, 0.7])
grid()
title("PA spectrum")
xlabel("288XXX GHz")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

figure()
NaCsPlot.plot_survival_data(datas[1][2], fmt="C0.-", label="Na+Cs/Na+Cs")
NaCsPlot.plot_survival_data(datas_cs[1][2], fmt="C1.-", label="Cs/Na+Cs")
NaCsPlot.plot_survival_data(datas_na[1][2], fmt="C2.-", label="Na/Na+Cs")
NaCsPlot.plot_survival_data(datas_na_na[1][2], fmt="C3.-", label="Na/Na")
legend(fontsize="small", ncol=2)
ylim([0.18, 0.9])
grid()
title("PA spectrum with \$\\sigma\$")
xlabel("288XXX GHz")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
