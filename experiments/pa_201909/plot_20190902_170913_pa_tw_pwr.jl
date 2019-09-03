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
    # datas = datas[[1:3:size(datas, 1) 2:3:size(datas, 1) 3:3:size(datas, 1)]]
    datas = datas[[1:2:size(datas, 1) 2:2:size(datas, 1)]]
    datas = NaCsData.map_params((i1, i2, v)->v[1], datas)
    return datas
end

const names_files = ["names_20190902_170913.mat"]
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
    ds[1] = NaCsData.map_params((i1, i2, v)->(v - 496) * 1000, ds[1])
end

function fit_survival(model, data, p0; plotx=nothing, plot_lo=nothing, plot_hi=nothing,
                      use_unc=true, plot_scale=1.1)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
    end
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        if plot_lo === nothing
            plot_lo = mid - span * plot_scale / 2
            if plot_lo * lo <= 0
                plot_lo = 0
            end
        end
        if plot_hi === nothing
            plot_hi = mid + span * plot_scale / 2
            if plot_hi * hi <= 0
                plot_hi = 0
            end
        end
        plotx = linspace(plot_lo, plot_hi, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios[:, 2], uncs[:, 2].^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios[:, 2], p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end

const prefix = joinpath(@__DIR__, "imgs", "data_20190902_170913_pa_tw_pwr")

# fit_5 = fit_survival(model_lorentzian, datas[1][:, 1], [0.7, 0.35, -20, 50])
# fit_10 = fit_survival(model_lorentzian, datas[1][:, 2], [0.7, 0.35, 850, 50])

figure()
NaCsPlot.plot_survival_data(datas[1][:, 1], fmt="C0.-", label="5 mW")
# plot(fit_5.plotx, fit_5.ploty, "C0")
NaCsPlot.plot_survival_data(datas[1][:, 2], fmt="C1.-", label="10 mW")
# plot(fit_10.plotx, fit_10.ploty, "C1")
# text(-175, 0.78, "\$\\Gamma_{5}=$(fit_5.uncs[4])\$ MHz", color="C0", fontsize="small")
# text(-175, 0.72, "\$f_{5}=$(fit_5.uncs[3])\$ MHz", color="C0", fontsize="small")
# text(750, 0.66, "\$\\Gamma_{10}=$(fit_10.uncs[4])\$ MHz", color="C1", fontsize="small")
# text(750, 0.6, "\$f_{10}=$(fit_10.uncs[3])\$ MHz", color="C1", fontsize="small")
grid()
ylim([0.27, 0.9])
# xlim([495.94, 496.16])
legend(fontsize="small", loc="upper right")
title("PA spectrum (NaCs/NaCs)")
xlabel("306496XXX MHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

pwrs = [5, 10]

figure(figsize=[12.6, 11.2])

for i in 1:2
    subplot(2, 2, i)
    NaCsPlot.plot_survival_data(datas_cs[1][:, i], fmt="C2.-", label="Cs/NaCs")
    NaCsPlot.plot_survival_data(datas_na[1][:, i], fmt="C1.-", label="Na/NaCs")
    NaCsPlot.plot_survival_data(datas[1][:, i], fmt="C0.-", label="NaCs/NaCs")
    grid()
    ylim([0.27, 0.9])
    # xlim([495.94, 496.16])
    legend(fontsize="small")
    title("PA spectrum (two-body $(pwrs[i]) mW)")
    xlabel("306496XXX MHz")
    ylabel("Two-body survival")
end

subplot(2, 2, 4)
NaCsPlot.plot_survival_data(datas_cscs[1][:, 1], fmt="C2.--", label="Cs 5 mW")
NaCsPlot.plot_survival_data(datas_cscs[1][:, 2], fmt="C2.:", label="Cs 10 mW")
NaCsPlot.plot_survival_data(datas_nana[1][:, 1], fmt="C1.--", label="Na 5 mW")
NaCsPlot.plot_survival_data(datas_nana[1][:, 2], fmt="C1.:", label="Na 10 mW")
grid()
ylim([0.27, 0.9])
# xlim([495.94, 496.16])
legend(fontsize="small", ncol=2)
title("PA spectrum (one-body)")
xlabel("306496XXX MHz")
ylabel("One-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_1b2b")

NaCsPlot.maybe_show()
