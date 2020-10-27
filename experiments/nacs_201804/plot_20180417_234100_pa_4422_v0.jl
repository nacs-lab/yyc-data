#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
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
            # sg = read(fd, "ScanGroup")
            # f = sg["base"]["params"]["fWavemeter1"]
            # ts = sg["scans"]["vars"]["params"]["TMergeWait"]
            f = pl[1]
            c = read(fd, "SingleAtomLogical") .!= 0
            return [f for i in 1:size(c, 3)], c
        end
        data = NaCsData.select_count(param, counts, selector)
        if !@isdefined datas
            datas = data
        else
            datas = [datas; data]
        end
    end
    # datas = datas[[1:2:size(datas, 1) 2:2:size(datas, 1)]]
    # datas = NaCsData.map_params((i1, i2, v)->v[1], datas)
    return datas
end

const names_files = ["names_20180417_234100.mat"
                     "names_20180418_071502.mat"]
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

for ds in [datas_cs, datas_na, datas, datas_cscs, datas_nana]
    ds[1] = NaCsData.map_params((i, v)->v * 11.576904 + 606.768489, ds[1])
    ds[2] = NaCsData.map_params((i, v)->v * 11.420860 + 608.195733, ds[2])
    # ds[1] = [ds[1]; ds[2]]
    # resize!(ds, 1)
end

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end

const prefix = joinpath(@__DIR__, "imgs", "data_20180417_234100_pa_4422_v0_j12")

# fit = fit_survival(model_lorentzian, datas[1], [0.8, .7, 492.35, 0.5])

figure()
axvline(693.33, color="C1", ls="--")
axvline(694.29, color="C1", ls="--")
axvline(698.03, color="C1", ls="--")
axvline(698.55, color="C1", ls="--")
axvline(705.39, color="C1", ls="--")
NaCsPlot.plot_survival_data(datas[1], fmt="C0.-")
NaCsPlot.plot_survival_data(datas[2], fmt="C0.-")
grid()
ylim([0, 0.95])
title("PA spectrum (two-body PA)")
xlabel("288XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
