#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200315_133410.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [([0, 0.06, 0.11, 0.16, 0.31, 0.46, 0.61] .* 1,
                [0, 0.06, 0.11, 0.16, 0.31, 0.46, 0.61] .* 1.5,
                [0, 0.06, 0.11, 0.16, 0.31, 0.46, 0.61] .* 3,
                [0, 0.06, 0.11, 0.16, 0.31, 0.46, 0.61] .* 8)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

fit_data(model, x, y, p0; kws...) =
    fit_data(model, x, y, nothing, p0; kws...)

function fit_data(model, params, ratios, uncs, p0;
                  plotx=nothing, plot_lo=nothing, plot_hi=nothing, plot_scale=1.1)
    use_unc = uncs !== nothing
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
        fit = curve_fit(model, params, ratios, uncs.^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios, p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

function fit_survival(model, data, p0; use_unc=true, kws...)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
        return fit_data(model, params, ratios[:, 2], uncs[:, 2], p0; kws...)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
        return fit_data(model, params, ratios[:, 2], p0; kws...)
    end
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const data_12 = datas_nacs[1][1]
const data_9 = datas_nacs[1][2]
const data_6 = datas_nacs[1][3]
const data_3 = datas_nacs[1][4]

const prefix = joinpath(@__DIR__, "imgs", "data_20200315_133410_raman_time_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
function model_expoff(x, p)
    p[1] .+ p[2] .* exp.(.-x ./ p[3])
end
# fit_t = fit_survival(model_expoff, data_t, [0.45, 0.25, 0.5])

figure()
NaCsPlot.plot_survival_data(data_12, fmt="C0.-", label="12 mW")
NaCsPlot.plot_survival_data(data_9, fmt="C1.-", label="9 mW")
NaCsPlot.plot_survival_data(data_6, fmt="C2.-", label="6 mW")
legend(fontsize="small")
title("288560 GHz")
xlim([0, 1.9])
grid()
xlabel("Raman time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
