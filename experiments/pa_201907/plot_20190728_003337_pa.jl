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

const names_files = ["names_20190728_003337.mat"]
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

const prefix = joinpath(@__DIR__, "imgs", "data_20190728_003337_pa")

fit_30 = fit_survival(model_lorentzian, datas[1][:, 1], [0.8, 0.3, 466.75, 0.05])
fit_75 = fit_survival(model_lorentzian, datas[1][:, 2], [0.8, 0.3, 466.75, 0.05])

figure()
NaCsPlot.plot_survival_data(datas[1][:, 1], fmt="C0.", label="0.3, 55 ms")
plot(fit_30.plotx, fit_30.ploty, "C0-")
NaCsPlot.plot_survival_data(datas[1][:, 2], fmt="C1.", label="0.75, 25 ms")
plot(fit_75.plotx, fit_75.ploty, "C1-")
text(466.54, 0.81, "\$\\Gamma_{0.3}=$(fit_30.uncs[4] * 1000)\$ MHz", color="C0")
text(466.54, 0.73, "\$\\Gamma_{0.75}=$(fit_75.uncs[4] * 1000)\$ MHz", color="C1")
grid()
ylim([0.4, 0.9])
legend(fontsize="small", loc="lower right")
title("PA spectrum (NaCs/NaCs)")
xlabel("307XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_2b")

figure()
NaCsPlot.plot_survival_data(datas_cs[1][:, 1], fmt="C2.-", label="Cs/NaCs")
NaCsPlot.plot_survival_data(datas_na[1][:, 1], fmt="C1.-", label="Na/NaCs")
NaCsPlot.plot_survival_data(datas[1][:, 1], fmt="C0.-", label="NaCs/NaCs")
grid()
ylim([0.4, 0.9])
legend(fontsize="small")
title("PA spectrum (two-body 0.30)")
xlabel("307XXX GHz")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_2b_30")

figure()
NaCsPlot.plot_survival_data(datas_cs[1][:, 2], fmt="C2.-", label="Cs/NaCs")
NaCsPlot.plot_survival_data(datas_na[1][:, 2], fmt="C1.-", label="Na/NaCs")
NaCsPlot.plot_survival_data(datas[1][:, 2], fmt="C0.-", label="NaCs/NaCs")
grid()
ylim([0.4, 0.9])
legend(fontsize="small")
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
ylim([0.4, 0.9])
legend(fontsize="small")
title("PA spectrum (one-body)")
xlabel("307XXX GHz")
ylabel("One-body survival")
NaCsPlot.maybe_save("$(prefix)_1b")

NaCsPlot.maybe_show()
