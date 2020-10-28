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
            f = sg["base"]["params"]["fWavemeter0"]
            amps = sg["base"]["vars"]["params"]["Merge"]["Cs"]["Wait"]["PTwzr"]
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
    datas = datas[[1:3:size(datas, 1) 2:3:size(datas, 1) 3:3:size(datas, 1)]]
    datas = NaCsData.map_params((i1, i2, v)->v[1], datas)
    return datas
end

const names_files = ["names_20190830_072607.mat"]
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

# for ds in [datas_cs, datas_na, datas, datas_cscs, datas_nana]
#     ds[1] = NaCsData.map_params((i1, i2, v)->v - 3000, ds[1])
#     ds[2] = NaCsData.map_params((i1, i2, v)->v - 3000, ds[2])
# end

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end

const prefix = joinpath(@__DIR__, "imgs", "data_20190830_072607_pa_tw_pwr")

fit_5 = fit_survival(model_lorentzian, datas[1][:, 2], [0.7, 0.35, 496.06, 0.1])

figure()
NaCsPlot.plot_survival_data(datas[1][:, 1], fmt="C0.-", label="1 mW")
NaCsPlot.plot_survival_data(datas[1][:, 2], fmt="C1.", label="5 mW")
plot(fit_5.plotx, fit_5.ploty, "C1-")
NaCsPlot.plot_survival_data(datas[1][:, 3], fmt="C2.-", label="10 mW")
text(495.945, 0.505, "\$\\Gamma_5=$(fit_5.uncs[4] * 1000)\$ MHz", color="C1")
text(495.945, 0.455, "\$f_5=$(fit_5.uncs[3])\$ GHz", color="C1")
grid()
ylim([0.45, 0.95])
xlim([495.94, 496.16])
legend(fontsize="small", loc="upper right")
title("PA spectrum (NaCs/NaCs)")
xlabel("306XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

pwrs = [1, 5, 10]

figure(figsize=[12.6, 11.2])

for i in 1:3
    subplot(2, 2, i)
    NaCsPlot.plot_survival_data(datas_cs[1][:, i], fmt="C2.-", label="Cs/NaCs")
    NaCsPlot.plot_survival_data(datas_na[1][:, i], fmt="C1.-", label="Na/NaCs")
    NaCsPlot.plot_survival_data(datas[1][:, i], fmt="C0.-", label="NaCs/NaCs")
    grid()
    ylim([0.45, 0.95])
    xlim([495.94, 496.16])
    legend(fontsize="small")
    title("PA spectrum (two-body $(pwrs[i]) mW)")
    xlabel("306XXX GHz")
    ylabel("Two-body survival")
end

subplot(2, 2, 4)
NaCsPlot.plot_survival_data(datas_cscs[1][:, 1], fmt="C2.-", label="Cs 1 mW")
NaCsPlot.plot_survival_data(datas_cscs[1][:, 2], fmt="C2.--", label="Cs 5 mW")
NaCsPlot.plot_survival_data(datas_cscs[1][:, 3], fmt="C2.:", label="Cs 10 mW")
NaCsPlot.plot_survival_data(datas_nana[1][:, 1], fmt="C1.-", label="Na 1 mW")
NaCsPlot.plot_survival_data(datas_nana[1][:, 2], fmt="C1.--", label="Na 5 mW")
NaCsPlot.plot_survival_data(datas_nana[1][:, 3], fmt="C1.:", label="Na 10 mW")
grid()
ylim([0.45, 0.95])
xlim([495.94, 496.16])
legend(fontsize="small", ncol=2)
title("PA spectrum (one-body)")
xlabel("306XXX GHz")
ylabel("One-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_1b2b")

NaCsPlot.maybe_show()
