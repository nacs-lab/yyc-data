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
            sg = read(fd, "ScanGroup")
            f = sg["base"]["params"]["fWavemeter1"]
            # amps = sg["base"]["vars"]["params"][2]["Merge"]["PA"]["Amp"]
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

const names_files = ["names_20200303_011858.mat",
                     "names_20200303_085714.mat"]
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
    # ds[1] = NaCsData.map_params((i, v)->(v - 306000), ds[1])
    ds[1] = [ds[1]; ds[2]]
    resize!(ds, 1)
end

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end

const prefix = joinpath(@__DIR__, "imgs", "data_20200303_011858_pa_3322")

# fit = fit_survival(model_lorentzian, datas[1], [0.8, 0.7, 492.35, 0.5])

figure()
NaCsPlot.plot_survival_data(datas_cs[1], fmt="C2.-", label="Cs/NaCs")
NaCsPlot.plot_survival_data(datas_na[1], fmt="C1.-", label="Na/NaCs")
NaCsPlot.plot_survival_data(datas[1], fmt="C0.-", label="NaCs/NaCs")
# plot(fit.plotx, fit.ploty, "C0-")
# text(490.82, 0.44, "\$\\Gamma=$(fit.uncs[4] * 1000)\$ MHz", color="C0", fontsize="small")
# text(490.82, 0.35, "\$f=$(fit.uncs[3])\$ GHz", color="C0", fontsize="small")
grid()
ylim([0, 1])
legend(fontsize="small")
title("PA spectrum (two-body PA)")
xlabel("288XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_2b")

figure()
NaCsPlot.plot_survival_data(datas_cscs[1], fmt="C2.-", label="Cs")
NaCsPlot.plot_survival_data(datas_nana[1], fmt="C1.-", label="Na")
grid()
ylim([0, 1])
# xlim([495.95, 612.29])
legend(fontsize="small")
title("PA spectrum (one-body)")
xlabel("288XXX GHz")
ylabel("One-body survival")
NaCsPlot.maybe_save("$(prefix)_1b")

NaCsPlot.maybe_show()
