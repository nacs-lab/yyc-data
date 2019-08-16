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
    datas = datas[[1:2:size(datas, 1) 2:2:size(datas, 1)]]
    datas = NaCsData.map_params((i1, i2, v)->v[1], datas)
    return datas
end

const names_files = ["names_20190815_234740.mat",
                     "names_20190816_081015.mat"]
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
    ds[1] = NaCsData.map_params((i1, i2, v)->v - 3000, ds[1])
    ds[2] = NaCsData.map_params((i1, i2, v)->v - 3000, ds[2])
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

const prefix = joinpath(@__DIR__, "imgs", "data_20190815_234740_pa")

fit_3 = fit_survival(model_lorentzian, [datas[1][:, 1]; datas[2][:, 1]], [0.8, 0.3, 255.45, 0.4])
fit_15 = fit_survival(model_lorentzian, [datas[1][:, 2]; datas[2][:, 2]], [0.8, 0.3, 255.45, 0.4])

figure()
NaCsPlot.plot_survival_data([datas[1][:, 1]; datas[2][:, 1]], fmt="C0.", label="3 mW")
plot(fit_3.plotx, fit_3.ploty, "C0-")
NaCsPlot.plot_survival_data([datas[1][:, 2]; datas[2][:, 2]], fmt="C1.", label="15 mW")
plot(fit_15.plotx, fit_15.ploty, "C1-")
text(255.0, 0.8, "\$\\Gamma_{3}=$(fit_3.uncs[4] * 1000)\$ MHz", color="C0")
text(255.3, 0.36, "\$\\Gamma_{15}=$(fit_15.uncs[4])\$ GHz", color="C1")
grid()
ylim([0.35, 0.95])
legend(fontsize="small", loc="upper right")
title("PA spectrum (NaCs/NaCs)")
xlabel("309XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_2b")

figure()
NaCsPlot.plot_survival_data([datas_cs[1][:, 1]; datas_cs[2][:, 1]], fmt="C2.-", label="Cs/NaCs")
NaCsPlot.plot_survival_data([datas_na[1][:, 1]; datas_na[2][:, 1]], fmt="C1.-", label="Na/NaCs")
NaCsPlot.plot_survival_data([datas[1][:, 1]; datas[2][:, 1]], fmt="C0.-", label="NaCs/NaCs")
grid()
ylim([0.35, 0.95])
legend(fontsize="small")
title("PA spectrum (two-body 3 mW)")
xlabel("309XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_2b_3")

figure()
NaCsPlot.plot_survival_data([datas_cs[1][:, 2]; datas_cs[2][:, 2]], fmt="C2.-", label="Cs/NaCs")
NaCsPlot.plot_survival_data([datas_na[1][:, 2]; datas_na[2][:, 2]], fmt="C1.-", label="Na/NaCs")
NaCsPlot.plot_survival_data([datas[1][:, 2]; datas[2][:, 2]], fmt="C0.-", label="NaCs/NaCs")
grid()
ylim([0.35, 0.95])
legend(fontsize="small")
title("PA spectrum (two-body 15 mW)")
xlabel("309XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_2b_15")

figure()
NaCsPlot.plot_survival_data([datas_cscs[1][:, 1]; datas_cscs[2][:, 1]], fmt="C2.-", label="Cs 3 mW")
NaCsPlot.plot_survival_data([datas_nana[1][:, 1]; datas_nana[2][:, 1]], fmt="C1.-", label="Na 3 mW")
NaCsPlot.plot_survival_data([datas_cscs[1][:, 2]; datas_cscs[2][:, 2]], fmt="C2.--", label="Cs 15 mW")
NaCsPlot.plot_survival_data([datas_nana[1][:, 2]; datas_nana[2][:, 2]], fmt="C1.--", label="Na 15 mW")
grid()
ylim([0.35, 0.95])
legend(fontsize="small")
title("PA spectrum (one-body)")
xlabel("309XXX GHz")
ylabel("One-body survival")
NaCsPlot.maybe_save("$(prefix)_1b")

NaCsPlot.maybe_show()
